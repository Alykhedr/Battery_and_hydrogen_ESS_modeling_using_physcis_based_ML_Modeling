% run_models.m — one switch to run B1 (SOC), B2 (+cycle), or B3 (+calendar)
clear; clc;

%% ---------------- Choose model here ----------------
% Options: 'B1' (energy only), 'B2' (cycle ageing), 'B3' (cycle+calendar)
model = 'B3';

% Optional: repeat the input year R times (e.g., long-horizon)
R = 1;

%% --------------- Load profile ----------------------
[time, Ppv, Pload] = load_profile_bess();
time  = time(:);  Ppv = Ppv(:);  Pload = Pload(:);

if R > 1
    % repeat time & signals seamlessly
    t1 = time; dt = t1(2)-t1(1);
    per = t1(end)-t1(1)+dt;
    offsets = (0:R-1)';
    time  = repmat(t1,R,1) + repelem(offsets*per, numel(t1));
    Ppv   = repmat(Ppv,  R, 1);
    Pload = repmat(Pload,R,1);
end

%% --------------- Params --------------
p = params(model);

%% --------------- Run the unified core --------------
out = core_loop(time, Ppv, Pload, p);

%% --------------- Plots -----------------------------
% Your existing plot: PV vs Load + SOC
plot_bess(out.time, Ppv, Pload, out.SOC, out.SOH);

% Extra plot for B2/B3: SOH vs EFC (kept inline like your old scripts)
if p.modes.do_cycle
    figure('Color','w');
    plot(out.EFC, 100*out.SOH, 'LineWidth',1.1); grid on;
    xlabel('EFC'); ylabel('SOH (%)');
    title(sprintf('SOH vs EFC — %s', p.meta.model));
end

%% --------------- Report ----------------------------
report(out);

%% --------------- Per-step contributions (no core edits) ---------------
if p.modes.do_cycle || p.modes.do_calendar
    % time axis for steps (start time of each interval)
    t_step = out.time(1:end-1);
    % per-step loss increments [percent-points]
    dSOH_cyc_pp = diff(out.loss.cycle);      % size N-1
    dSOH_cal_pp = diff(out.loss.calendar);   % size N-1

    % reconstruct achieved currents from terminal power & efficiencies
    % P_bess(k) = + (V*Idis/1000)*eta_dis  -  (V*Ich/1000)/eta_ch
    P = out.P_bess(:);                       % kW, +dis / -chg
    V = p.V_nom;                             % V
    eta_dis = p.eta_dis;  eta_ch = p.eta_ch;

    Idis_A = zeros(size(P));                 % A, achieved discharge current
    Ich_A  = zeros(size(P));                 % A, achieved charge current
    dis_mask = (P >= 0);
    chg_mask = ~dis_mask;
    Idis_A(dis_mask) = (1000*P(dis_mask)) ./ (V*eta_dis);
    Ich_A(chg_mask)  = (-1000*P(chg_mask)) .* eta_ch ./ V;

    % achieved C-rates
    C_dch = Idis_A / p.Q_nom;
    C_ch  = Ich_A  / p.Q_nom;

    % choose which EFC to display:
    EFC_start = out.EFC(1:end-1);    % EFC at start of step k
    EFC_end   = out.EFC(2:end);      % EFC after the step (k+1)

    % temperature column (constant here; adapt if you later make T(t))
    T_degC = p.T_degC * ones(size(P));

    % build tidy table
    T = table( ...
        t_step,         ... % interval start time
        P,              ... % kW
        Idis_A, Ich_A,  ... % A
        C_dch,  C_ch,   ... % C-rate
        T_degC,         ...
        EFC_start, EFC_end, ...
        dSOH_cyc_pp, dSOH_cal_pp, ...
        'VariableNames', {'t_start','P_bess_kW','I_dis_A','I_ch_A','C_dch','C_ch','T_degC','EFC_start','EFC_end','dSOH_cyc_pp','dSOH_cal_pp'});


    outdir = 'results';
    if ~exist(outdir,'dir'), mkdir(outdir); end
    writetable(T, fullfile(outdir, sprintf('step_contrib_%s.csv', p.meta.model)));
end
%% ----- Save per-step CSV robustly (absolute path + diagnostics) -----
% Resolve this script's folder (works for .m scripts)
thisFile = mfilename('fullpath');
if isempty(thisFile)
    % Fallback to current folder if running from Live Script/Command Window
    baseDir = pwd;
else
    baseDir = fileparts(thisFile);
end

outdir = fullfile(baseDir, 'results');

% Make folder (with status + message)
[mkOK, mkMsg, mkMsgID] = mkdir(outdir);
if ~mkOK
    warning('Could not create results folder: %s (%s). Using pwd instead.', mkMsg, mkMsgID);
    outdir = fullfile(pwd, 'results');
    [mkOK2, mkMsg2] = mkdir(outdir);
    if ~mkOK2
        error('Failed to create fallback results folder: %s', mkMsg2);
    end
end

% Timestamped filename to avoid overwrite / caching confusions
ts      = datestr(now, 'yyyymmdd_HHMMSS');
csvName = sprintf('step_contrib_%s_%s.csv', p.meta.model, ts);
csvPath = fullfile(outdir, csvName);

% Try writing, report any error explicitly
try
    writetable(T, csvPath);
catch ME
    fprintf(2, 'writetable failed: %s\n', ME.message);
    rethrow(ME);
end

% Verify existence, and print where it went
if exist(csvPath, 'file') == 2
    fprintf('✓ Saved per-step CSV: %s\n', csvPath);
else
    warning('CSV not found after write. Check write permissions for: %s', outdir);
end
