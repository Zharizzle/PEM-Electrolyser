%% Automatic calibration for Simulink -- reads user's Excel file, picks
%% the best model on its own
%
% WHAT THIS DOES
%   1. Reads the user's Excel file (current density/voltage data).
%   2. Reads the membrane thickness straight off the "Membrane Thickness"
%      Constant block on the model -- a genuine user input on the model
%      itself, typed in exactly like n_cells or Active Area, never edited
%      here.
%   3. Tries every "hold one temperature out, check the rest" test,
%      the same way the separate diagnostic script does, but silently --
%      just to decide which works better: Linear or Arrhenius.
%   4. Once it knows which form wins, it re-fits ONE final time using
%      ALL the temperatures (not holding anything back), since the
%      hold-out testing was only needed to judge which form to trust,
%      not to throw away real data from the final answer.
%   5. Writes the final numbers into the workspace, where Simulink
%      picks them up automatically.
%
% HOW TO CONNECT THIS TO SIMULINK (one-time setup):
%   Model Settings (gear icon) -> Callbacks -> InitFcn -> paste:
%       run_model_calibration
%   This makes Simulink run this script automatically every time you
%   press Run -- so it always uses whatever Excel file is currently
%   sitting in the folder, and whatever membrane thickness is currently
%   typed into the model, with no separate manual step.
%
% WHERE TO PUT THE EXCEL FILE:
%   Same folder as this script, named exactly: user_data.xlsx
%   Columns required, exact names: Temperature_C,
%   CurrentDensity_Acm2, Voltage_V (membrane thickness comes from the
%   Membrane Thickness block on the model, not the spreadsheet)

clear; clc;

R = 8.314;
F = 96485.3;
a_water = 3;
data_file = 'user_data.xlsx';

if ~isfile(data_file)
    error('Could not find %s in the current folder. Place your data file there and try again.', data_file);
end

Tbl = readtable(data_file);
all_T_C = unique(Tbl.Temperature_C)';

%% ================= Membrane thickness: read from the model ================= %%
% Membrane thickness is a user input on the model, not in this script.
% Type the value (in micrometres) directly into the "Membrane Thickness"
% Constant block on the canvas, exactly like n_cells or Active Area. This
% script just reads whatever number is currently sitting in that block,
% so there is nowhere else to change it and nothing that can fall out of
% sync between the calibration and the live model.
mem_block = [bdroot, '/Membrane Thickness'];   % <-- EDIT the block name here if yours differs
t_mem_um = str2double(get_param(mem_block, 'Value'));
if isnan(t_mem_um) || t_mem_um <= 0
    error('Could not read a valid membrane thickness from %s. Check the block value.', mem_block);
end
%% ============================================================================ %%
t_mem = t_mem_um * 1e-4;   % converts to cm for internal use

if length(all_T_C) < 4
    fprintf('WARNING: only %d temperature(s) in this data. With fewer than 4,\n', length(all_T_C));
    fprintf('         the correlation cannot be properly tested -- results may be unreliable.\n');
end

fprintf('Found %d temperature(s) in %s: %s\n', length(all_T_C), data_file, mat2str(all_T_C));

%% ================= i_lim bounds: scale to the data, not a fixed number ================= %%
% i_lim must be at least twice the highest current density actually
% tested (Section 3.4) -- a fitted i_lim below the max tested current is
% non-physical, since current cannot be measured above the limiting
% current. This mirrors the bound used in the batch calibration-fit
% script exactly, so the two scripts can never disagree on this.
max_i_tested = max(Tbl.CurrentDensity_Acm2);
i_lim_margin = 2.0;
lb = [0.05, 0.05, 1e-6, 1e-6, max(0.3, i_lim_margin*max_i_tested)];
ub = [2.0,  2.0,  0.5,  0.5,  max(5.0, lb(5)+1.0)];
x0 = [0.5, 0.5, 1e-4, 1e-3, lb(5) + 0.5];
%% ========================================================================================= %%
to_bounded   = @(z) lb + (ub - lb) ./ (1 + exp(-z));
to_unbounded = @(x) -log((ub - lb)./(x - lb) - 1);

%% --- Step 1: silently test Linear vs Arrhenius via leave-one-out ---
lin_scores = []; arr_scores = [];
if length(all_T_C) >= 4
    for h = 1:length(all_T_C)
        holdout = all_T_C(h);
        calib_T = all_T_C(all_T_C ~= holdout);
        [lin_R2, arr_R2] = test_one_holdout(Tbl, calib_T, holdout, t_mem, R, F, a_water, lb, ub, x0, to_bounded, to_unbounded);
        lin_scores(end+1) = lin_R2; %#ok<AGROW>
        arr_scores(end+1) = arr_R2; %#ok<AGROW>
    end
    fprintf('Mean held-out R^2 -- Linear: %.4f | Arrhenius: %.4f\n', mean(lin_scores), mean(arr_scores));
    if mean(arr_scores) >= mean(lin_scores)
        use_form = 'Arrhenius';
    else
        use_form = 'Linear';
    end
else
    fprintf('Not enough temperatures to test forms -- defaulting to Linear.\n');
    use_form = 'Linear';
end
fprintf('Selected form: %s\n', use_form);

%% --- Step 2: final fit using ALL temperatures (no holding back) ---
T_K = all_T_C + 273.15;
n_T = length(T_K);
expData = cell(1, n_T);
for k = 1:n_T
    rows = abs(Tbl.Temperature_C - all_T_C(k)) < 0.5;
    expData{k} = [Tbl.CurrentDensity_Acm2(rows), Tbl.Voltage_V(rows)];
end

fitted = zeros(n_T, 5);
for k = 1:n_T
    Tk = T_K(k);
    i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
    objective = @(z) sum((V_exp - arrayfun(@(ii) ...
        PEMVoltage_fit(ii, Tk, to_bounded_elem(z,lb,ub,1), to_bounded_elem(z,lb,ub,2), ...
            to_bounded_elem(z,lb,ub,3), to_bounded_elem(z,lb,ub,4), to_bounded_elem(z,lb,ub,5), ...
            R, F, t_mem, a_water), i_exp)).^2);
    z_fit = fminsearch(objective, to_unbounded(x0), optimset('Display','off','MaxFunEvals',3000));
    fitted(k,:) = to_bounded(z_fit);
end
i_lim_final = mean(fitted(:,5));
T_ref_final = median(T_K);

lin_coefs = zeros(4,2); arr_coefs = zeros(4,2);
for p = 1:4
    lin_coefs(p,:) = polyfit(T_K, fitted(:,p)', 1);
    y = fitted(:,p)';
    obj = @(x) sum((y - x(1)*exp(-x(2)/R*(1./T_K - 1/T_ref_final))).^2);
    arr_coefs(p,:) = fminsearch(obj, [mean(y), 5000], optimset('Display','off','MaxFunEvals',3000));
end

%% --- Push to workspace for Simulink ---
if strcmp(use_form, 'Linear')
    form_type = 0; %#ok<NASGU>
    kinetic_coefs = reshape(lin_coefs', 1, []); %#ok<NASGU>
else
    form_type = 1; %#ok<NASGU>
    kinetic_coefs = reshape(arr_coefs', 1, []); %#ok<NASGU>
end
T_ref = T_ref_final; %#ok<NASGU>
i_lim = i_lim_final; %#ok<NASGU>

%% --- R^2 for the selected form, always printed ---
% With >=4 temperatures this is the proper held-out R^2 already printed
% above. With <4 temperatures the held-out loop never ran (2 calibration
% points can't meaningfully test a 2-parameter correlation -- Section
% 3.5), so instead evaluate the selected correlation back at the
% temperatures it was built from. This is calibration-fit R^2, not
% held-out, and is labelled as such so it is never read as a real test.
if length(all_T_C) < 4
    fit_R2v = zeros(n_T, 1);
    for k = 1:n_T
        Tk = T_K(k); i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
        if strcmp(use_form, 'Linear')
            aa=polyval(lin_coefs(1,:),Tk); ac=polyval(lin_coefs(2,:),Tk);
            j0a=polyval(lin_coefs(3,:),Tk); j0c=polyval(lin_coefs(4,:),Tk);
        else
            aa=arr_coefs(1,1)*exp(-arr_coefs(1,2)/R*(1/Tk-1/T_ref_final));
            ac=arr_coefs(2,1)*exp(-arr_coefs(2,2)/R*(1/Tk-1/T_ref_final));
            j0a=arr_coefs(3,1)*exp(-arr_coefs(3,2)/R*(1/Tk-1/T_ref_final));
            j0c=arr_coefs(4,1)*exp(-arr_coefs(4,2)/R*(1/Tk-1/T_ref_final));
        end
        aa=min(max(aa,lb(1)),ub(1)); ac=min(max(ac,lb(2)),ub(2)); j0a=min(max(j0a,lb(3)),ub(3)); j0c=min(max(j0c,lb(4)),ub(4));
        Vp = arrayfun(@(ii) PEMVoltage_fit(ii,Tk,aa,ac,j0a,j0c,i_lim_final,R,F,t_mem,a_water), i_exp);
        fit_R2v(k) = 1 - sum((V_exp-Vp).^2)/sum((V_exp-mean(V_exp)).^2);
    end
    fprintf('Calibration-fit R^2 (NOT held-out -- see note below) for %s: %.4f\n', use_form, mean(fit_R2v));
end

%% --- Print the selected correlation, always (even when <4 temperatures) ---
% The held-out test only runs with >=4 temperatures; below that, use_form
% defaults to Linear without ever being tested, and the coefficients were
% previously never shown. Print them regardless, flagged as unreliable
% when the held-out test could not run, so nothing is silently hidden.
param_names = {'alpha_a','alpha_c','j0_a','j0_c'};

if strcmp(use_form, 'Linear')
    for p = 1:4
        fprintf('  %s = %.4e*T + %.4e;\n', param_names{p}, lin_coefs(p,1), lin_coefs(p,2));
    end
else
    fprintf('  T_ref = %.2f;\n', T_ref_final);
    for p = 1:4
        fprintf('  %s = %.4e * exp(-%.4e/R * (1/T - 1/T_ref));\n', param_names{p}, arr_coefs(p,1), arr_coefs(p,2));
    end
end
fprintf('  %% i_lim: %.4f\n', i_lim_final);
if length(all_T_C) < 4
    fprintf('  NOTE: fewer than 4 temperatures -- this correlation could not be held-out tested and may be unreliable.\n');
    fprintf('  NOTE: the R^2 above is calibration-fit (built on these same points), not held-out -- see Section 3.5.\n');
end

fprintf('\nDone. Simulink will now use: %s form, i_lim = %.4f', use_form, i_lim);

%% --- Local functions ---
function val = to_bounded_elem(z, lb, ub, idx)
    val = lb(idx) + (ub(idx) - lb(idx)) / (1 + exp(-z(idx)));
end

function [lin_R2, arr_R2] = test_one_holdout(Tbl, calib_T, holdout_T, t_mem, R, F, a_water, lb, ub, x0, to_bounded, to_unbounded)
    T_K = calib_T + 273.15;
    n_T = length(T_K);
    expData = cell(1, n_T);
    for k = 1:n_T
        rows = abs(Tbl.Temperature_C - calib_T(k)) < 0.5;
        expData{k} = [Tbl.CurrentDensity_Acm2(rows), Tbl.Voltage_V(rows)];
    end
    fitted = zeros(n_T, 5);
    for k = 1:n_T
        Tk = T_K(k);
        i_exp = expData{k}(:,1); V_exp = expData{k}(:,2);
        objective = @(z) sum((V_exp - arrayfun(@(ii) ...
            PEMVoltage_fit(ii, Tk, to_bounded_elem(z,lb,ub,1), to_bounded_elem(z,lb,ub,2), ...
                to_bounded_elem(z,lb,ub,3), to_bounded_elem(z,lb,ub,4), to_bounded_elem(z,lb,ub,5), ...
                R, F, t_mem, a_water), i_exp)).^2);
        z_fit = fminsearch(objective, to_unbounded(x0), optimset('Display','off','MaxFunEvals',2000));
        fitted(k,:) = to_bounded(z_fit);
    end
    i_lim_h = mean(fitted(:,5));
    T_ref_h = median(T_K);
    lin_coefs = zeros(4,2); arr_coefs = zeros(4,2);
    for p = 1:4
        lin_coefs(p,:) = polyfit(T_K, fitted(:,p)', 1);
        y = fitted(:,p)';
        obj = @(x) sum((y - x(1)*exp(-x(2)/R*(1./T_K - 1/T_ref_h))).^2);
        arr_coefs(p,:) = fminsearch(obj, [mean(y), 5000], optimset('Display','off','MaxFunEvals',2000));
    end
    Tk_h = holdout_T + 273.15;
    rows = abs(Tbl.Temperature_C - holdout_T) < 0.5;
    i_hold = Tbl.CurrentDensity_Acm2(rows); V_hold = Tbl.Voltage_V(rows);

    aa=polyval(lin_coefs(1,:),Tk_h); ac=polyval(lin_coefs(2,:),Tk_h);
    j0a=polyval(lin_coefs(3,:),Tk_h); j0c=polyval(lin_coefs(4,:),Tk_h);
    aa=min(max(aa,lb(1)),ub(1)); ac=min(max(ac,lb(2)),ub(2)); j0a=min(max(j0a,lb(3)),ub(3)); j0c=min(max(j0c,lb(4)),ub(4));
    Vp = arrayfun(@(ii) PEMVoltage_fit(ii,Tk_h,aa,ac,j0a,j0c,i_lim_h,R,F,t_mem,a_water), i_hold);
    lin_R2 = 1 - sum((V_hold-Vp).^2)/sum((V_hold-mean(V_hold)).^2);

    aa=arr_coefs(1,1)*exp(-arr_coefs(1,2)/R*(1/Tk_h-1/T_ref_h));
    ac=arr_coefs(2,1)*exp(-arr_coefs(2,2)/R*(1/Tk_h-1/T_ref_h));
    j0a=arr_coefs(3,1)*exp(-arr_coefs(3,2)/R*(1/Tk_h-1/T_ref_h));
    j0c=arr_coefs(4,1)*exp(-arr_coefs(4,2)/R*(1/Tk_h-1/T_ref_h));
    aa=min(max(aa,lb(1)),ub(1)); ac=min(max(ac,lb(2)),ub(2)); j0a=min(max(j0a,lb(3)),ub(3)); j0c=min(max(j0c,lb(4)),ub(4));
    Vp = arrayfun(@(ii) PEMVoltage_fit(ii,Tk_h,aa,ac,j0a,j0c,i_lim_h,R,F,t_mem,a_water), i_hold);
    arr_R2 = 1 - sum((V_hold-Vp).^2)/sum((V_hold-mean(V_hold)).^2);
end