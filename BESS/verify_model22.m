% verify_B2_closedform_naumann.m
% Cycle-aging (B2, Montes LFP) verification without any power profile.
% We evaluate SOH(FEC) in closed form and compare to Naumann cycle datasets.
%
% Needs in path: params.m (with B2 constants)
% Input: one Naumann MAT file with X_Axis_Data_Mat, Y_Axis_Data_Mat, Legend_Vec

clear; clc;

%% ---- choose a MAT file ----
fn = 'xDOD_1C1C_40_C_Capacity_CC_CV_FEC.mat';   % e.g., 'xCyC_80DOD_40_C_Capacity_CC_CV_FEC.mat'

S = load(fn);
% try common variable names; error if missing
assert(isfield(S,'X_Axis_Data_Mat') && isfield(S,'Y_Axis_Data_Mat') && isfield(S,'Legend_Vec'), ...
    'Expected X_Axis_Data_Mat, Y_Axis_Data_Mat, Legend_Vec in %s', fn);

FEC_mat = S.X_Axis_Data_Mat;      % columns = series
Cap_mat = S.Y_Axis_Data_Mat;      % capacity (often normalized to 1 at FEC=0)
Leg     = S.Legend_Vec;           % cellstr

fprintf('=== B2 (closed-form) verification — %s ===\n', fn);

%% ---- parameters (B2) ----
p = params('B2');                 % cycle ON, calendar OFF (not used here)
a  = p.a_montes;

%% ---- loop series ----
rows = {};
figure('Color','w'); hold on; grid on;
for i = 1:numel(Leg)
    label = strtrim(Leg{i});
    % skip CC+CV (CV dwell not in B2)
    if contains(lower(label),'cc+cv') || contains(lower(label),'cc_cv')
        fprintf('  [skip] %s (CC+CV)\n', label);
        continue;
    end

    F = FEC_mat(:,i);  Q = Cap_mat(:,i);
    m = ~isnan(F) & ~isnan(Q);
    F = F(m); Q = Q(m);

    if isempty(F), continue; end

    % measured SOH (normalize to first)
    SOH_meas = Q / Q(1);

    % parse conditions from legend/filename
    [Tdeg, soc_ctr, dod_pct, Cch, Cdch] = parse_legend(label, fn);
    mSOC = soc_ctr/100;

    % Montes multipliers (DOD in percent, T in K)
    T_K = Tdeg + 273.15;
    fT  = exp( p.kT * (T_K - p.Tref_K) / T_K );
    fD  = exp( p.kDODc * dod_pct );                 % NOTE: DOD in percent

    delta = p.kcyc * fT * fD * exp(p.kCch*Cch) * exp(p.kCdch*Cdch) ...
          * (1 + p.kmSOC * mSOC * ((1 - mSOC)/(2*p.mSOCref)));

    % closed-form SOH(F) = 1 - (delta/100) * F^a
    SOH_model = 1 - (delta/100) * (F.^a);
    if isfield(p,'deg') && isfield(p.deg,'SOH_min')
        SOH_model = max(SOH_model, p.deg.SOH_min);
    end

    % metrics
    rmse = sqrt(mean((SOH_model - SOH_meas).^2,'omitnan'));

    i80_d = find(SOH_meas <= 0.80, 1, 'first');
    i80_m = find(SOH_model<= 0.80, 1, 'first');
    EFC80_data  = iff(isempty(i80_d), NaN, F(i80_d));
    EFC80_model = iff(isempty(i80_m), NaN, F(i80_m));

    % print + stash
    fprintf('  %-55s  T=%2.0f°C  DOD=%3.0f%%  Cch=%.2gC  Cdch=%.2gC  RMSE=%.3f  EFC80 data=%s  model=%s\n', ...
        clip(label,55), Tdeg, dod_pct, Cch, Cdch, rmse, numfmt(EFC80_data), numfmt(EFC80_model));
    rows(end+1,:) = {label, Tdeg, dod_pct, Cch, Cdch, rmse, EFC80_data, EFC80_model};

    % overlay
    plot(F, 100*SOH_meas, 'o', 'DisplayName',[label ' (data)']);
    plot(F, 100*SOH_model, '-', 'DisplayName',[label ' (B2)']);
end

xlabel('FEC'); ylabel('SOH (%)');
title(sprintf('Naumann CC series vs B2 (closed form) — %s', fn), 'Interpreter','none');
legend('Location','bestoutside');

if ~isempty(rows)
    Tsum = cell2table(rows, 'VariableNames', ...
        {'Series','T_degC','DOD_pct','Cch','Cdch','RMSE','EFC80_data','EFC80_model'});
    disp('Summary:'); disp(Tsum);
else
    fprintf('No CC series processed.\n');
end

%% --------- helpers (local functions) ----------
function [Tdeg, soc_pct, dod_pct, Cch, Cdch] = parse_legend(label, fallback_fn)
    Tdeg    = getnum(label, '(\d+)\s*°C', NaN);
    soc_pct = getnum(label, '(\d+)\s*%SOC', NaN);
    dod_pct = getnum(label, '(\d+)\s*%DOD', NaN);
    cr = regexp(label, '(\d+(?:\.\d+)?)C[_\-xX](\d+(?:\.\d+)?)C', 'tokens', 'once');
    if isempty(cr)
        cr = regexp(label, '(\d+(?:\.\d+)?)C.*?(\d+(?:\.\d+)?)C', 'tokens', 'once');
    end
    if isempty(cr), Cch = NaN; Cdch = NaN;
    else, Cch = str2double(cr{1}); Cdch = str2double(cr{2});
    end
    % fallbacks
    if isnan(Tdeg),    Tdeg    = getnum(fallback_fn, '(\d+)\s*_?C', 25); end
    if isnan(soc_pct), soc_pct = 50; end
    if isnan(dod_pct), dod_pct = getnum(fallback_fn, '(\d+)\s*DOD', 100); end
    if isnan(Cch),     Cch = 1.0; end
    if isnan(Cdch),    Cdch = 1.0; end
end

function v = getnum(str, pat, default)
    t = regexp(str, pat, 'tokens', 'once');
    if isempty(t), v = default; else, v = str2double(t{1}); end
end

function s = numfmt(x)
    if isnan(x), s = 'NaN'; else, s = sprintf('%.0f', x); end
end

function y = iff(c,a,b), if c, y=a; else, y=b; end, end

function s = clip(txt,n)
    if numel(txt)<=n, s=txt; else, s=[txt(1:n-1) '…']; end
end
