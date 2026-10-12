%% Batch Stage 1 / Stage 2 calibration on ALL temperatures, for ALL datasets

clear; clc;

R = 8.314;
F = 96485.3;
a_water = 3;
data_file = 'all_polarisation_data.xlsx';

Tbl_all = readtable(data_file);
sources_available = unique(Tbl_all.Source, 'stable');

activation_region = [0.01, 0.10];
ohmic_region      = [0.20, 0.60];

results = struct();

for s = 1:numel(sources_available)
    setup_source = sources_available{s};
    setup_rows = strcmp(Tbl_all.Source, setup_source);
    t_mem = Tbl_all.Thickness_um(find(setup_rows, 1)) * 1e-4;   % um -> cm
    all_T_C = unique(Tbl_all.Temperature_C(setup_rows))';       % row vector
    calib_T_C = all_T_C;
    calib_T_K = calib_T_C + 273.15;
    n_T = length(calib_T_C);

    fprintf('\n========== %s (%.0f um) ==========\n', strrep(setup_source,'_',' '), t_mem*1e4);
    fprintf('Calibration temperatures (ALL, no holdout): %s C\n', mat2str(calib_T_C));

    expData = cell(1, n_T);
    for k = 1:n_T
        rows = setup_rows & abs(Tbl_all.Temperature_C - calib_T_C(k)) < 0.5;
        expData{k} = [Tbl_all.CurrentDensity_Acm2(rows), Tbl_all.Voltage_V(rows)];
    end

    max_i_tested = max(cellfun(@(d) max(d(:,1)), expData));
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
        z_fit = fminsearch(objective, z0, optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000));
        x_fit = to_bounded(z_fit);
        fitted(k,:) = x_fit;
        V_pred = arrayfun(@(ii) PEMVoltage_fit(ii, Tk, x_fit(1),x_fit(2),x_fit(3),x_fit(4),x_fit(5), R,F,t_mem,a_water), i_exp);
        SS_res = sum((V_exp-V_pred).^2); SS_tot = sum((V_exp-mean(V_exp)).^2);
        fit_R2(k) = 1 - SS_res/SS_tot;
        fprintf('  %d C: alpha_a=%.4f alpha_c=%.4f j0_a=%.2e j0_c=%.2e i_lim=%.3f (R^2=%.4f)\n', ...
            calib_T_C(k), x_fit(1),x_fit(2),x_fit(3),x_fit(4),x_fit(5), fit_R2(k));
    end
    i_lim_final = mean(fitted(:,5));
    fprintf('  i_lim_final (mean across all temperatures) = %.4f A/cm^2\n', i_lim_final);

    param_names = {'alpha_a','alpha_c','j0_a','j0_c'};
    lin_coefs = zeros(4,2); arr_coefs = zeros(4,2);
    T_ref = median(calib_T_K);
    for p = 1:4
        lin_coefs(p,:) = polyfit(calib_T_K, fitted(:,p)', 1);
        y = fitted(:,p)';
        obj = @(x) sum((y - x(1)*exp(-x(2)/R*(1./calib_T_K - 1/T_ref))).^2);
        arr_coefs(p,:) = fminsearch(obj, [mean(y), 5000], optimset('Display','off','MaxFunEvals',5000));
    end

    lin_R2v = zeros(n_T,1); arr_R2v = zeros(n_T,1);
    for k = 1:n_T
        Tk = calib_T_K(k); i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
        aa=polyval(lin_coefs(1,:),Tk); ac=polyval(lin_coefs(2,:),Tk);
        j0a=polyval(lin_coefs(3,:),Tk); j0c=polyval(lin_coefs(4,:),Tk);
        Vp = arrayfun(@(ii) PEMVoltage_fit(ii,Tk,aa,ac,j0a,j0c,i_lim_final,R,F,t_mem,a_water), i_exp);
        lin_R2v(k) = 1 - sum((V_exp-Vp).^2)/sum((V_exp-mean(V_exp)).^2);

        aa=arr_coefs(1,1)*exp(-arr_coefs(1,2)/R*(1/Tk-1/T_ref));
        ac=arr_coefs(2,1)*exp(-arr_coefs(2,2)/R*(1/Tk-1/T_ref));
        j0a=arr_coefs(3,1)*exp(-arr_coefs(3,2)/R*(1/Tk-1/T_ref));
        j0c=arr_coefs(4,1)*exp(-arr_coefs(4,2)/R*(1/Tk-1/T_ref));
        Vp = arrayfun(@(ii) PEMVoltage_fit(ii,Tk,aa,ac,j0a,j0c,i_lim_final,R,F,t_mem,a_water), i_exp);
        arr_R2v(k) = 1 - sum((V_exp-Vp).^2)/sum((V_exp-mean(V_exp)).^2);
    end

    fprintf('  T(C)   Linear   Arrhenius\n');
    for k = 1:n_T
        fprintf('  %4d   %.4f   %.4f\n', calib_T_C(k), lin_R2v(k), arr_R2v(k));
    end
    fprintf('  Mean: linear=%.4f | Arrhenius=%.4f\n', mean(lin_R2v), mean(arr_R2v));

    if mean(arr_R2v) >= mean(lin_R2v)
        use_form = 'Arrhenius';
    else
        use_form = 'Linear';
    end
    fprintf('  --> Auto-selected form (by calibration-fit R^2): %s\n', use_form);

    fprintf('  --- Copy-paste-ready lines for PEMVoltage.m ---\n');
    if strcmp(use_form, 'Linear')
        for p = 1:4
            fprintf('  %s = %.4e*T + %.4e;\n', param_names{p}, lin_coefs(p,1), lin_coefs(p,2));
        end
    else
        fprintf('  T_ref = %.2f;\n', T_ref);
        for p = 1:4
            fprintf('  %s = %.4e * exp(-%.4e/R * (1/T - 1/T_ref));\n', param_names{p}, arr_coefs(p,1), arr_coefs(p,2));
        end
    end
    fprintf('  %% i_lim: %.4f\n', i_lim_final);

    % --- Plot: per-temperature calibration fits, all temperatures ---
    figure('Visible','off'); hold on; grid on;
    colors_calib = jet(n_T);
    for k = 1:n_T
        Tk = calib_T_K(k);
        i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
        i_sweep_c = linspace(0.001, max(i_exp)*1.02, 100);
        V_sweep_c = arrayfun(@(ii) PEMVoltage_fit(ii, Tk, fitted(k,1), fitted(k,2), fitted(k,3), fitted(k,4), fitted(k,5), R, F, t_mem, a_water), i_sweep_c);
        plot(i_sweep_c, V_sweep_c, 'Color', colors_calib(k,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%d C (R^2=%.3f)', calib_T_C(k), fit_R2(k)));
        plot(i_exp, V_exp, 'o', 'Color', colors_calib(k,:), 'MarkerFaceColor', colors_calib(k,:), ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    xlabel('Current Density (A/cm^2)'); ylabel('Cell Voltage (V)');
    title(sprintf('%s: Per-Temperature Calibration Fits (All Temperatures)', strrep(setup_source,'_',' ')));
    legend('show', 'Location', 'northwest');
    fname = sprintf('fig4_1_%s_calib.svg', matlab.lang.makeValidName(setup_source));
    exportgraphics(gcf, fname, 'ContentType', 'vector');
    fprintf('  Figure saved: %s\n', fname);
    close(gcf);

    % Results
    key = matlab.lang.makeValidName(setup_source);
    results.(key).T_C = calib_T_C;
    results.(key).fitted = fitted;
    results.(key).fit_R2 = fit_R2;
    results.(key).i_lim_final = i_lim_final;
    results.(key).use_form = use_form;
    results.(key).lin_coefs = lin_coefs;
    results.(key).arr_coefs = arr_coefs;
    results.(key).T_ref = T_ref;
end

fprintf('\nDone. All figures saved to the current folder. Results struct available as ''results'' in the workspace.\n');

function val = to_bounded_elem(z, lb, ub, idx)
    val = lb(idx) + (ub(idx) - lb(idx)) / (1 + exp(-z(idx)));
end