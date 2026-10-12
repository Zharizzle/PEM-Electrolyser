%% Block-cross validation sweep, per dataset

clear; clc;
apply_clamp = true;

R = 8.314;
F = 96485.3;
n_cells = 1;
a_water = 3;
data_file = 'all_polarisation_data.xlsx';

Tbl = readtable(data_file);
sources_available = unique(Tbl.Source, 'stable');

summary = struct('source',{},'holdout_T',{},'form',{},'held_out_R2',{},'held_out_RMSE',{},'clamp_flag',{},'clamped_params',{});
ds_summary = struct('source',{},'form',{},'mean_lin',{},'mean_arr',{},'n_T',{});

for s = 1:length(sources_available)
    setup_source = sources_available{s};
    setup_rows = strcmp(Tbl.Source, setup_source);
    t_mem = Tbl.Thickness_um(find(setup_rows,1)) * 1e-4;
    all_T_C = unique(Tbl.Temperature_C(setup_rows))';   % row vector

    n_holdouts = length(all_T_C);
    fprintf('\n========== %s (%d um) ==========\n', strrep(setup_source,'_',' '), t_mem*1e4);

    clear H;
    figure('Name', strrep(setup_source,'_',' '), 'Position', [50 50 1100 250*n_holdouts]);

    for h = 1:n_holdouts
        holdout_T_C = all_T_C(h);
        calib_T_C = all_T_C(all_T_C ~= holdout_T_C);
        calib_T_K = calib_T_C + 273.15;
        n_T = length(calib_T_C);
        T_ref = median(calib_T_K);

        expData = cell(1, n_T);
        for k = 1:n_T
            rows = setup_rows & abs(Tbl.Temperature_C - calib_T_C(k)) < 0.5;
            expData{k} = [Tbl.CurrentDensity_Acm2(rows), Tbl.Voltage_V(rows)];
        end

 
        max_i_tested = max(Tbl.CurrentDensity_Acm2(setup_rows));
        i_lim_margin = 2.0;
        lb = [0.05, 0.05, 1e-6, 1e-6, max(0.3, i_lim_margin*max_i_tested)];
        ub = [2.0,  2.0,  0.5,  0.5,  max(5.0, lb(5)+1.0)];
        x0 = [0.5, 0.5, 1e-4, 1e-3, lb(5) + 0.5];
        to_bounded   = @(z) lb + (ub - lb) ./ (1 + exp(-z));
        to_unbounded = @(x) -log((ub - lb)./(x - lb) - 1);
        z0 = to_unbounded(x0);

        fitted = zeros(n_T, 5);
        fit_R2 = zeros(n_T, 1);
        for k = 1:n_T
            Tk = calib_T_K(k);
            i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
            objective = @(z) sum((V_exp - arrayfun(@(ii) ...
                PEMVoltage_fit(ii, Tk, to_bounded_elem(z,lb,ub,1), to_bounded_elem(z,lb,ub,2), ...
                    to_bounded_elem(z,lb,ub,3), to_bounded_elem(z,lb,ub,4), to_bounded_elem(z,lb,ub,5), ...
                    R, F, t_mem, a_water), i_exp)).^2);
            z_fit = fminsearch(objective, z0, optimset('Display','off','MaxFunEvals',2000));  % same limits as run_model_calibration.m (test_one_holdout)
            x_fit = to_bounded(z_fit);
            fitted(k,:) = x_fit;
            V_pred = arrayfun(@(ii) PEMVoltage_fit(ii, Tk, x_fit(1),x_fit(2),x_fit(3),x_fit(4),x_fit(5), R,F,t_mem,a_water), i_exp);
            SS_res = sum((V_exp-V_pred).^2); SS_tot = sum((V_exp-mean(V_exp)).^2);
            fit_R2(k) = 1 - SS_res/SS_tot;
        end
        i_lim_final = mean(fitted(:,5));

        lin_coefs = zeros(4,2); arr_coefs = zeros(4,2);
        for p = 1:4
            lin_coefs(p,:) = polyfit(calib_T_K, fitted(:,p)', 1);
            y = fitted(:,p)';
            obj = @(x) sum((y - x(1)*exp(-x(2)/R*(1./calib_T_K - 1/T_ref))).^2);
            arr_coefs(p,:) = fminsearch(obj, [mean(y), 5000], optimset('Display','off','MaxFunEvals',2000));
        end

        Tk_hold = holdout_T_C + 273.15;
        hold_rows = setup_rows & abs(Tbl.Temperature_C - holdout_T_C) < 0.5;
        i_hold = Tbl.CurrentDensity_Acm2(hold_rows);
        V_hold = Tbl.Voltage_V(hold_rows);

        aa_lin=polyval(lin_coefs(1,:),Tk_hold); ac_lin=polyval(lin_coefs(2,:),Tk_hold);
        j0a_lin=polyval(lin_coefs(3,:),Tk_hold); j0c_lin=polyval(lin_coefs(4,:),Tk_hold);
        clamped_lin = {}; clamped_arr = {};
        pnames = {'alpha_a','alpha_c','j0_a','j0_c'};
        vals = [aa_lin,ac_lin,j0a_lin,j0c_lin];
        for pp=1:4
            if vals(pp)<=lb(pp) || vals(pp)>=ub(pp)
                clamped_lin{end+1} = [pnames{pp} '(lin)']; %#ok<AGROW>
            end
        end
        if apply_clamp
            aa_lin = min(max(aa_lin, lb(1)), ub(1)); ac_lin = min(max(ac_lin, lb(2)), ub(2));
            j0a_lin = min(max(j0a_lin, lb(3)), ub(3)); j0c_lin = min(max(j0c_lin, lb(4)), ub(4));
        end
        V_hold_pred_lin = arrayfun(@(ii) PEMVoltage_fit(ii, Tk_hold, aa_lin, ac_lin, j0a_lin, j0c_lin, i_lim_final, R, F, t_mem, a_water), i_hold);
        R2_hold_lin = 1 - sum((V_hold-V_hold_pred_lin).^2)/sum((V_hold-mean(V_hold)).^2);

        aa_arr=arr_coefs(1,1)*exp(-arr_coefs(1,2)/R*(1/Tk_hold-1/T_ref));
        ac_arr=arr_coefs(2,1)*exp(-arr_coefs(2,2)/R*(1/Tk_hold-1/T_ref));
        j0a_arr=arr_coefs(3,1)*exp(-arr_coefs(3,2)/R*(1/Tk_hold-1/T_ref));
        j0c_arr=arr_coefs(4,1)*exp(-arr_coefs(4,2)/R*(1/Tk_hold-1/T_ref));
        vals = [aa_arr,ac_arr,j0a_arr,j0c_arr];
        for pp=1:4
            if vals(pp)<=lb(pp) || vals(pp)>=ub(pp)
                clamped_arr{end+1} = [pnames{pp} '(arr)']; %#ok<AGROW>
            end
        end
        if apply_clamp
            aa_arr = min(max(aa_arr, lb(1)), ub(1)); ac_arr = min(max(ac_arr, lb(2)), ub(2));
            j0a_arr = min(max(j0a_arr, lb(3)), ub(3)); j0c_arr = min(max(j0c_arr, lb(4)), ub(4));
        end
        V_hold_pred_arr = arrayfun(@(ii) PEMVoltage_fit(ii, Tk_hold, aa_arr, ac_arr, j0a_arr, j0c_arr, i_lim_final, R, F, t_mem, a_water), i_hold);
        R2_hold_arr = 1 - sum((V_hold-V_hold_pred_arr).^2)/sum((V_hold-mean(V_hold)).^2);

        % store everything pass 2 needs for this holdout
        H(h).holdout_T_C = holdout_T_C; H(h).calib_T_C = calib_T_C; H(h).calib_T_K = calib_T_K;
        H(h).n_T = n_T; H(h).expData = expData; H(h).fitted = fitted; H(h).i_lim_final = i_lim_final;
        H(h).Tk_hold = Tk_hold; H(h).i_hold = i_hold; H(h).V_hold = V_hold;
        H(h).lin = [aa_lin, ac_lin, j0a_lin, j0c_lin]; H(h).arr = [aa_arr, ac_arr, j0a_arr, j0c_arr];
        H(h).Vp_lin = V_hold_pred_lin; H(h).Vp_arr = V_hold_pred_arr;
        H(h).R2_lin = R2_hold_lin; H(h).R2_arr = R2_hold_arr;
        H(h).clamped_lin = clamped_lin; H(h).clamped_arr = clamped_arr;
    end

    %% --- Choose the correlation form ONCE for this dataset ---
    mean_R2_lin = mean([H.R2_lin]);
    mean_R2_arr = mean([H.R2_arr]);
    if n_holdouts >= 4
        if mean_R2_arr >= mean_R2_lin
            use_form = 'Arrhenius';
        else
            use_form = 'Linear';
        end
        fprintf('  Mean held-out R^2 -- Linear: %.4f | Arrhenius: %.4f\n', mean_R2_lin, mean_R2_arr);
    else
        use_form = 'Linear';
        fprintf('  Only %d temperatures (<4): form not testable, defaulting to Linear.\n', n_holdouts);
        fprintf('  (For reference, mean held-out R^2 -- Linear: %.4f | Arrhenius: %.4f)\n', mean_R2_lin, mean_R2_arr);
    end
    fprintf('  --> Selected form for %s: %s (used for every holdout below)\n', strrep(setup_source,'_',' '), use_form);
    ds_summary(end+1) = struct('source',setup_source,'form',use_form,'mean_lin',mean_R2_lin,'mean_arr',mean_R2_arr, ...
        'n_T',n_holdouts); %#ok<SAGROW>

    %% --- Pass 2: report each holdout with the selected form ---
    for h = 1:n_holdouts
        holdout_T_C = H(h).holdout_T_C; calib_T_C = H(h).calib_T_C; calib_T_K = H(h).calib_T_K;
        n_T = H(h).n_T; expData = H(h).expData; fitted = H(h).fitted; i_lim_final = H(h).i_lim_final;
        Tk_hold = H(h).Tk_hold; i_hold = H(h).i_hold; V_hold = H(h).V_hold;
        if strcmp(use_form, 'Arrhenius')
            p4 = H(h).arr; V_hold_pred = H(h).Vp_arr; clamped_params = H(h).clamped_arr;
        else
            p4 = H(h).lin; V_hold_pred = H(h).Vp_lin; clamped_params = H(h).clamped_lin;
        end
        aa = p4(1); ac = p4(2); j0a = p4(3); j0c = p4(4);
        clamp_hit = ~isempty(clamped_params);
        fprintf('  [%d C holdout] held-out R^2: Linear=%.4f, Arrhenius=%.4f (reporting %s)\n', ...
            holdout_T_C, H(h).R2_lin, H(h).R2_arr, use_form);
        overall_R2 = 1 - sum((V_hold-V_hold_pred).^2)/sum((V_hold-mean(V_hold)).^2);
        overall_RMSE = sqrt(mean((V_hold-V_hold_pred).^2));

        %% --- Theory-derived region boundaries: sweep the model's own eta
        %% breakdown to find where each term actually dominates, instead of
        %% assuming a fixed literature window applies to every case. ---
        i_sweep_reg = linspace(0.001, max(i_hold), 300);
        eta_act_sweep = zeros(size(i_sweep_reg));
        eta_ohm_sweep = zeros(size(i_sweep_reg));
        eta_diff_sweep = zeros(size(i_sweep_reg));
        for si = 1:length(i_sweep_reg)
            [eta_act_sweep(si), eta_ohm_sweep(si), eta_diff_sweep(si)] = ...
                compute_etas(i_sweep_reg(si), Tk_hold, aa, ac, j0a, j0c, i_lim_final, R, F, t_mem, a_water);
        end
        [~, dominant] = max([eta_act_sweep; eta_ohm_sweep; eta_diff_sweep], [], 1);
        % dominant: 1=activation, 2=ohmic, 3=diffusion, per sweep point

        act_idx = i_sweep_reg(dominant==1);
        ohm_idx = i_sweep_reg(dominant==2);
        diff_idx = i_sweep_reg(dominant==3);

        if ~isempty(act_idx)
            activation_region = [min(act_idx), max(act_idx)];
        else
            activation_region = [NaN, NaN];
        end
        ohmic_region_is_fallback = false;
        if ~isempty(ohm_idx)
            ohmic_region = [min(ohm_idx), max(ohm_idx)];
        else
            % Degenerate case: activation never stops dominating within the
            % tested range -- almost always because a clamped alpha (at
            % extrapolation) artificially inflates eta_act everywhere.
            % Rather than report nothing, fall back to the top half of the
            % tested range, clearly flagged as a fallback rather than a
            % genuine theory-derived crossover.
            ohmic_region = [0.5*max(i_hold), max(i_hold)];
            ohmic_region_is_fallback = true;
        end
        if ~isempty(diff_idx)
            diffusion_region_start = min(diff_idx);
        else
            diffusion_region_start = Inf;
        end
        % True diffusion dominance requires i close to i_lim, but i_lim was
        % deliberately set >= 2x the tested current range (to avoid the
        % non-physical sharp elbow artifact), so that crossover is never
        % reached by real data -- diffusion_region_start above will
        % essentially always be Inf. Rather than report NaN for every
        % single case, this is redefined as the top 20% of the TESTED
        % current density range -- an honest "high-current behaviour"
        % check, not a claim of genuine mass-transport-limited dominance.
        high_current_region_start = 0.8 * max(i_hold);
        diffusion_region_start = high_current_region_start;

        fprintf('  Theory-derived regions: activation=[%.3f-%.3f], ohmic=[%.3f-%.3f]\n', ...
            activation_region(1), activation_region(2), ohmic_region(1), ohmic_region(2));
        fprintf('  Diffusion dominance not reached within tested range (i_lim >= 2x max tested, by design).\n');
        fprintf('  Reporting high-current region (top 20%% of tested range, starts %.3f A/cm^2) instead.\n', diffusion_region_start);

        % --- Activation region: Tafel slope comparison ---
        act_mask = i_hold >= activation_region(1) & i_hold <= activation_region(2);
        fallback_used_act = false;
        if sum(act_mask) < 2 && length(i_hold) >= 3
            % Sparse experimental sampling can miss a narrow theory-derived
            % window entirely even though the region genuinely exists --
            % fall back to the 3 held-out points nearest the window,
            % rather than reporting nothing.
            [~, sort_idx] = sort(abs(i_hold - mean(activation_region)));
            act_mask = false(size(i_hold));
            act_mask(sort_idx(1:min(3,length(i_hold)))) = true;
            fallback_used_act = true;
        end
        tafel_note = 'insufficient points';
        tafel_diff_pct = NaN;
        if sum(act_mask) >= 2
            i_act = i_hold(act_mask); V_act_exp = V_hold(act_mask); V_act_pred = V_hold_pred(act_mask);
            p_exp_act = polyfit(log(i_act), V_act_exp, 1);
            p_pred_act = polyfit(log(i_act), V_act_pred, 1);
            tafel_diff_abs = (p_pred_act(1) - p_exp_act(1)) * log(10);   % V/decade
            if abs(p_exp_act(1)) > 1e-4
                tafel_diff_pct = 100*(p_pred_act(1) - p_exp_act(1))/p_exp_act(1);
                pctstr = sprintf('%.1f%% diff', tafel_diff_pct);
            else
                tafel_diff_pct = NaN;
                pctstr = 'pct undefined, ref slope ~0';
            end
            fbnote = ''; if fallback_used_act, fbnote = ' [FALLBACK: nearest points, window had <2]'; end
            tafel_note = sprintf('exp slope=%.4f V/decade, model=%.4f V/decade (abs diff=%.4f V/decade, %s)%s', ...
                p_exp_act(1)*log(10), p_pred_act(1)*log(10), tafel_diff_abs, pctstr, fbnote);
        end
        fprintf('    Activation region [%.4f-%.4f A/cm^2]: %d held-out points -- %s\n', activation_region(1), activation_region(2), sum(act_mask), tafel_note);

        % --- Ohmic region: linear slope (ASR) comparison ---
        % IMPORTANT: comparing the model's *exact* isolated ohmic term
        % (t_mem/sigma alone) against a *local linear fit* of the real
        % total voltage curve is not a fair like-for-like comparison --
        % the real curve's local slope still carries a small residual
        % contribution from eta_act (which does not vanish once ohmic
        % dominates, only shrinks) and possibly eta_diff, so it is not
        % measuring pure ohmic behaviour the way the exact value does.
        % This was found to explain a consistent, physically-sensible
        % pattern of mismatch even in datasets with excellent overall R^2
        % (e.g. Kumar). The fix: compare the model's OWN local-fit slope
        % (fit to its smooth, noise-free predicted curve over the same
        % window) against the experimental local-fit slope -- both sides
        % now measuring the same thing. The exact physics value is still
        % reported separately, as a useful secondary check of how much
        % residual curvature the local window is picking up.
        a_c_ohm = min(max(a_water, 0), 3);
        if a_c_ohm <= 1
            lambda_ohm = 0.043 + 17.81*a_c_ohm - 39.85*a_c_ohm^2 + 36.0*a_c_ohm^3;
        else
            lambda_ohm = 14 + 1.4*(a_c_ohm - 1);
        end
        sigma_ohm = max((0.005139*lambda_ohm - 0.00326) * exp(1268*(1/303 - 1/Tk_hold)), 1e-6);
        model_ASR_exact = t_mem / sigma_ohm;   % pure ohmic term alone, for reference only

        ohm_mask = i_hold >= ohmic_region(1) & i_hold <= ohmic_region(2);
        ohmic_note = 'insufficient points for a slope estimate (need >=2)';
        ohmic_diff_pct = NaN;
        ohmic_diff_abs = NaN;
        ohmic_R2 = NaN;
        if sum(ohm_mask) >= 2
            i_ohm = i_hold(ohm_mask); V_ohm_exp = V_hold(ohm_mask); V_ohm_pred = V_hold_pred(ohm_mask);
            p_exp_ohm = polyfit(i_ohm, V_ohm_exp, 1);
            p_pred_ohm = polyfit(i_ohm, V_ohm_pred, 1);   % model's LOCAL slope -- now the primary comparison
            ohmic_diff_abs = p_pred_ohm(1) - p_exp_ohm(1);   % Ohm.cm^2, both sides now local fits, fair comparison
            if abs(p_exp_ohm(1)) > 1e-3
                ohmic_diff_pct = 100*ohmic_diff_abs/p_exp_ohm(1);
                pctstr2 = sprintf('%.1f%% diff', ohmic_diff_pct);
            else
                ohmic_diff_pct = NaN;
                pctstr2 = 'pct undefined, ref slope ~0';
            end
            % Note: intercept/offset comparison removed from reporting
            % per request -- ASR (slope) comparison remains the primary
            % ohmic diagnostic.
            fbflag = ''; if ohmic_region_is_fallback, fbflag = ' [FALLBACK REGION]'; end
            exact_gap = model_ASR_exact - p_pred_ohm(1);
            ohmic_note = sprintf('exp ASR (local fit)=%.4f Ohm.cm^2, model ASR (local fit)=%.4f Ohm.cm^2 (abs diff=%.4f, %s) | pure ohmic term=%.4f (residual curvature in window=%.4f)%s', ...
                p_exp_ohm(1), p_pred_ohm(1), ohmic_diff_abs, pctstr2, model_ASR_exact, exact_gap, fbflag);
        elseif sum(ohm_mask) == 1
            % A single point can't define a slope at all, but we can still
            % report how close the model's prediction is to that one point.
            i_ohm = i_hold(ohm_mask); V_ohm_exp = V_hold(ohm_mask); V_ohm_pred = V_hold_pred(ohm_mask);
            ohmic_R2 = -abs(V_ohm_pred - V_ohm_exp);   % not a real R^2, single point -- report as a signed gap instead
            ohmic_note = sprintf('only 1 held-out point in region -- no slope possible; point-wise voltage gap=%.4fV', V_ohm_pred - V_ohm_exp);
        end
        fbtag = ''; if ohmic_region_is_fallback, fbtag = ' [FALLBACK]'; end
        fprintf('    Ohmic region [%.4f-%.4f A/cm^2]%s: %d held-out points -- %s\n', ohmic_region(1), ohmic_region(2), ...
            fbtag, sum(ohm_mask), ohmic_note);

        % --- Diffusion region: RMSE (slope not meaningful here) ---
        diff_mask = i_hold >= diffusion_region_start;
        diff_note = 'no held-out points reach this region';
        diff_RMSE = NaN;
        if sum(diff_mask) >= 1
            diff_RMSE = sqrt(mean((V_hold(diff_mask)-V_hold_pred(diff_mask)).^2));
            diff_note = sprintf('RMSE=%.4f V (n=%d points)', diff_RMSE, sum(diff_mask));
        end
        fprintf('    Diffusion region [starts %.4f A/cm^2, max held-out i=%.4f]: %d points -- %s\n', diffusion_region_start, max(i_hold), sum(diff_mask), diff_note);


        summary(end+1) = struct('source',setup_source,'holdout_T',holdout_T_C,'form',use_form, ...
            'held_out_R2',overall_R2,'held_out_RMSE',overall_RMSE,'clamp_flag',clamp_hit,'clamped_params',{clamped_params}); %#ok<SAGROW>

        clampnote = ''; if clamp_hit, clampnote = sprintf(' [CLAMP: %s]', strjoin(unique(clamped_params), ', ')); end
        fprintf('  Holdout %d C: %s form, held-out R^2=%.4f, RMSE=%.4f%s\n', ...
            holdout_T_C, use_form, overall_R2, overall_RMSE, clampnote);

        % --- Left subplot: calibration fits (this holdout's calibration set) ---
        subplot(n_holdouts, 2, 2*h-1); hold on; grid on;
        colors_c = jet(n_T);
        for k = 1:n_T
            Tk = calib_T_K(k);
            i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
            i_sweep_c = linspace(0.001, max(i_exp)*1.02, 80);
            V_sweep_c = arrayfun(@(ii) PEMVoltage_fit(ii, Tk, fitted(k,1),fitted(k,2),fitted(k,3),fitted(k,4),fitted(k,5), R,F,t_mem,a_water), i_sweep_c);
            plot(i_sweep_c, V_sweep_c, 'Color', colors_c(k,:), 'LineWidth', 1);
            plot(i_exp, V_exp, 'o', 'Color', colors_c(k,:), 'MarkerSize', 3, 'HandleVisibility','off');
        end
        title(sprintf('Block-Cross Calibration (Holdout = %d C): T = %s', holdout_T_C, mat2str(calib_T_C)), 'FontSize', 8);
        xlabel('A/cm^2','FontSize',7); ylabel('V','FontSize',7);

        % --- Right subplot: held-out validation, with theory-derived regions shaded ---
        subplot(n_holdouts, 2, 2*h); hold on; grid on;
        i_sweep = linspace(0.001, max(i_hold)*1.02, 150);
        V_sweep = arrayfun(@(ii) PEMVoltage_fit(ii, Tk_hold, aa, ac, j0a, j0c, i_lim_final, R, F, t_mem, a_water), i_sweep);

        yl_lo = min([V_hold; V_sweep']) * 0.98;
        yl_hi = max([V_hold; V_sweep']) * 1.02;

        % Regions are still computed from theory (eta_act/eta_ohm/eta_diff
        % crossover, above) -- no colour shading drawn on the plot anymore,
        % per updated preference. The region boundaries and their
        % diagnostics are reported as text only (see regiontxt below and
        % console output above).
        if ohmic_region_is_fallback
            fprintf('  NOTE: theory-derived ohmic region was degenerate (activation dominates the whole\n');
            fprintf('        tested range -- likely a clamped alpha). Using top-half-of-range fallback instead.\n');
        end

        plot(i_sweep, V_sweep, 'b-', 'LineWidth', 1.3);
        plot(i_hold, V_hold, 'ro', 'MarkerSize', 4);
        ylim([yl_lo yl_hi]);

        clamptxt = ''; if clamp_hit, clamptxt = sprintf(' [%s]', strjoin(unique(clamped_params), ',')); end
        title(sprintf('Block-Cross Validation, Held-Out %d C (%s, R^2 = %.3f)', holdout_T_C, use_form, overall_R2), 'FontSize', 8);
        xlabel('A/cm^2','FontSize',7); ylabel('V','FontSize',7);

        % Region-specific error metrics, printed directly on the plot
        actstr = 'n/a'; if isfinite(tafel_diff_pct), actstr = sprintf('%.0f%%', tafel_diff_pct); elseif exist('tafel_diff_abs','var') && isfinite(tafel_diff_abs), actstr = sprintf('%.3fV/dec', tafel_diff_abs); end
        ohmstr = 'n/a'; if isfinite(ohmic_diff_pct), ohmstr = sprintf('%.0f%%', ohmic_diff_pct); elseif exist('ohmic_diff_abs','var') && isfinite(ohmic_diff_abs), ohmstr = sprintf('%.4fOhm.cm^2', ohmic_diff_abs); elseif exist('ohmic_R2','var') && isfinite(ohmic_R2), ohmstr = sprintf('1pt gap=%.3fV', -ohmic_R2); end
        regiontxt = sprintf('Act: %s | Ohm: %s | HighCurrent RMSE: %.3fV', actstr, ohmstr, diff_RMSE);
        xl = xlim;
        text(xl(1)+0.02*(xl(2)-xl(1)), yl_lo+0.05*(yl_hi-yl_lo), regiontxt, 'FontSize', 6, 'Color', [0.3 0.3 0.3]);
    end
    sgtitle(sprintf('%s: Block-Cross Validation Sweep', strrep(setup_source,'_',' ')));
end

%% --- Overall summary table ---
fprintf('\n\n=== SUMMARY: held-out R^2 by dataset and holdout temperature ===\n');
fprintf('%-22s %8s %10s %10s %10s %-20s\n', 'Source', 'Holdout', 'Form', 'R^2', 'RMSE', 'Clamped params');
for i = 1:length(summary)
    if summary(i).clamp_flag
        cp_str = strjoin(unique(summary(i).clamped_params), ',');
    else
        cp_str = '-';
    end
    fprintf('%-22s %8d %10s %10.4f %10.4f %-20s\n', strrep(summary(i).source,'_',' '), ...
        summary(i).holdout_T, summary(i).form, summary(i).held_out_R2, summary(i).held_out_RMSE, cp_str);
end

fprintf('\n=== Selected form per dataset (compare with Table 4.1 mean held-out R^2) ===\n');
fprintf('%-22s %5s %10s %12s %12s\n', 'Source', 'nT', 'Form', 'mean R2 Lin', 'mean R2 Arr');
for i = 1:length(ds_summary)
    fprintf('%-22s %5d %10s %12.4f %12.4f\n', strrep(ds_summary(i).source,'_',' '), ds_summary(i).n_T, ...
        ds_summary(i).form, ds_summary(i).mean_lin, ds_summary(i).mean_arr);
end

function val = to_bounded_elem(z, lb, ub, idx)
    val = lb(idx) + (ub(idx) - lb(idx)) / (1 + exp(-z(idx)));
end

function [eta_act, eta_ohm, eta_diff] = compute_etas(I, T, alpha_a, alpha_c, j0_a, j0_c, i_lim, R, F, t_mem, a)
%COMPUTE_ETAS  Same physics as PEMVoltage_fit.m, but returns the three
%   overpotential terms separately instead of the summed voltage --
%   used here purely to find where each term dominates, for
%   theory-derived region boundaries (not a duplicate of the main model,
%   since PEMVoltage_fit.m's single-output signature is relied on
%   elsewhere via arrayfun and should not be changed).
    a_c = min(max(a, 0), 3);
    if a_c <= 1
        lambda = 0.043 + 17.81*a_c - 39.85*a_c^2 + 36.0*a_c^3;
    else
        lambda = 14 + 1.4*(a_c - 1);
    end
    sigma_raw = (0.005139*lambda - 0.00326) * exp(1268*(1/303 - 1/T));
    sigma = max(sigma_raw, 1e-6);

    eta_act = 0; eta_ohm = 0; eta_diff = 0;
    if I > 0
        eta_act = (R*T/(alpha_a*F)) * asinh(I/(2*j0_a)) + (R*T/(alpha_c*F)) * asinh(I/(2*j0_c));
        eta_ohm = I * (t_mem/sigma);
        i_safe = min(I, 0.999*i_lim);
        eta_diff = (R*T/F) * log(i_lim/(i_lim - i_safe));
    end
end