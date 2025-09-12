% verify_B3_calendar.m — Validate calendar aging (no loads) vs Naumann data
% Uses DischargeCapacity.xlsx (Naumann 2020 JPS calendar dataset).
% Assumes your codebase has: params('B3'), core_loop, calendar_aging.

clear; clc;

%% --- CONFIG: choose a few representative test points by header name
% Examples from the sheet: 'TP_25°C,50%SOC', 'TP_40°C,100%SOC', 'TP_40°C,0%SOC', 'TP_60°C,50%SOC'
TP_LIST = { 'TP_0°C,50%SOC','TP_25°C,0%SOC','TP_40°C,0%SOC','TP_25°C,50%SOC', 'TP_40°C,100%SOC', 'TP_10°C,50%SOC', 'TP_40°C,50%SOC' };

% Path to the calendar capacity file
XLS_FILE = 'DischargeCapacity.xlsx';   % place the Naumann file in the working dir

%% --- Load table (robust to minor header variations)
T = readtable(XLS_FILE, 'VariableNamingRule','preserve');  
% Find the time column (hours)
time_col = find(contains(string(T.Properties.VariableNames), "Storage", 'IgnoreCase',true) ...
              & contains(string(T.Properties.VariableNames), "hour",   'IgnoreCase',true), 1);
if isempty(time_col)
    error('Could not find "Storage time / hours" column in %s', XLS_FILE);
end
t_h = T{:, time_col};
t_h = t_h(:);

%% --- Loop over requested test points
fprintf('=== B3 Calendar Verification (Naumann calendar dataset) ===\n');
for k = 1:numel(TP_LIST)
    tp_name = TP_LIST{k};
    % Try exact header; if not found, try forgiving match (digits °C and %SOC)
    col_idx = find(strcmp(T.Properties.VariableNames, tp_name), 1);
    if isempty(col_idx)
        col_idx = find_smart(T.Properties.VariableNames, tp_name);
    end
    if isempty(col_idx)
        warning('Test point "%s" not found. Skipping.', tp_name); %#ok<WNTAG>
        continue;
    end

    Q = T{:, col_idx};           % discharge capacity in Ah
    if all(isnan(Q)), warning('All NaN in "%s". Skipping.', tp_name); continue; end

    % Build measured SOH from first non-NaN
    i0  = find(~isnan(Q), 1, 'first');
    Q0  = Q(i0);
    SOH_meas = Q ./ Q0;

    % Parse T (°C) and SOC setpoint from the header label
    [T_degC, SOC_set] = parse_tp_label(tp_name);   % SOC_set in [0..1]

    % --- Run model with calendar ONLY, zero loads, fixed SOC
    % Time vector = dataset timestamps; zero PV and zero load keep SOC constant
    time = hours(t_h - t_h(1));
    Ppv   = zeros(numel(time),1);
    Pload = zeros(numel(time),1);

    p = params('B3');                 % calendar ON by default
    p.modes.do_cycle    = false;      % ensure cycle OFF
    p.modes.do_calendar = true;
    p.T_degC = T_degC;

    % Keep the SOC fixed at the setpoint by giving a narrow window around it
    margin = 1e-6;                    % practically fixed
    p.SOC_min = max(0, SOC_set - margin);
    p.SOC_max = min(1, SOC_set + margin);
    p.SOC0    = SOC_set;

    out = core_loop(time, Ppv, Pload, p);

    % Align model SOH to measured samples (same time grid already)
    SOH_pred = out.SOH(:);

    % Compute metrics on overlapping valid points
    mask = ~isnan(SOH_meas) & ~isnan(SOH_pred);
    rmse = sqrt(mean((SOH_pred(mask) - SOH_meas(mask)).^2));
    fin_err = SOH_pred(find(mask,1,'last')) - SOH_meas(find(mask,1,'last'));

    fprintf('%-22s  T=%2.0f°C  SOC=%3.0f%%  RMSE=%.003f  FinalΔ=%.003f\n', ...
        tp_name, T_degC, round(100*SOC_set), rmse, fin_err);

    % Optional quick visual (comment out if you want print-only)
    figure('Color','w','Name',tp_name); 
    plot(t_h, 100*SOH_meas,'o', 'DisplayName','Measured'); hold on; grid on;
    plot(t_h, 100*SOH_pred,'-', 'DisplayName','Model');
    xlabel('Storage time (h)'); ylabel('SOH (%)'); title(tp_name); legend('Location','southwest');
end

%% ----------------- helpers -----------------

function idx = find_smart(varnames, label)
% Try to locate a column by matching "<num>°C,<num>%SOC" inside varnames
    idx = [];
    [Tc, SOCc] = parse_tp_label(label);
    patC = sprintf('%d°C', round(Tc));
    patS = sprintf('%d%%SOC', round(100*SOCc));
    for i = 1:numel(varnames)
        vn = varnames{i};
        if contains(vn, matlab.lang.makeValidName(patC,'ReplacementStyle','delete')) ...
        && contains(vn, matlab.lang.makeValidName(patS,'ReplacementStyle','delete'))
            idx = i; return;
        end
    end
end

function [Tc, SOCf] = parse_tp_label(s)
% Parse 'TP_40°C,100%SOC' → Tc=40, SOCf=1.00
    % strip non-ASCII degree if needed
    s = strrep(s,'°','°'); 
    % find temperature
    cidx = regexp(s, '(\d+)\s*°C', 'tokens', 'once');
    if isempty(cidx), cidx = regexp(s, '(\d+)\s*C', 'tokens', 'once'); end
    Tc = str2double(cidx{1});
    % find SOC percent
    pidx = regexp(s, '(\d+)\s*%SOC', 'tokens', 'once');
    SOCf = str2double(pidx{1})/100;
end
