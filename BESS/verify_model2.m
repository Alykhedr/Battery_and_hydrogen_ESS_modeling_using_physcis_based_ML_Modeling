% verify_B2_sweeps_simple.m
% Montes LFP (B2) verification sweeps with clean prints (no plots).
% - Full-window DOD (0..100%), calendar OFF
% - CC charge then CC discharge each cycle
% - EFC reported at SOH = 80% (or "censored" if not reached)

clear; clc;

% --- Params (cycle ON, calendar OFF) ---
p = params('B2');

% Convenience
P1C_kW = (p.V_nom * p.Q_nom)/1000;   % 1C pack power
fprintf('=== B2 Verification (Montes LFP) ===\n');
fprintf('Pack: V_nom=%.1f V, Q_nom=%.0f Ah  →  1C power = %.1f kW\n', ...
        p.V_nom, p.Q_nom, P1C_kW);

% Helper to compute EFC@80% SOH for given (T, Cch, Cdch)
efc80 = @(TT, Cch, Cdch) efc_at_80pct(p, TT, Cch, Cdch);

%% 1) Temperature sweep at 1C/1C (expect life ↓ as T ↑)
temps = [30 40 50];
life_T = nan(size(temps));
for i = 1:numel(temps)
    life_T(i) = efc80(temps(i), 1.0, 1.0);
end
fprintf('\n-- Temperature sweep (Cch=1C, Cdch=1C) --\n');
for i = 1:numel(temps)
    if isnan(life_T(i))
        fprintf('T=%2d°C : EFC@80%% > limit (censored)\n', temps(i));
    else
        fprintf('T=%2d°C : EFC@80%% = %.0f\n', temps(i), life_T(i));
    end
end
if all(isfinite(life_T))
    trendT = all(diff(life_T) <= 0);  % non-increasing vs T
    fprintf('Check (life should decrease with T): %s\n', tern(trendT,'OK','WARN'));
end

%% 2) Discharge C-rate sweep at 35°C (expect strongest effect)
Cdchs = [2 3 4 5 6];
life_Cd = nan(size(Cdchs));
for i = 1:numel(Cdchs)
    life_Cd(i) = efc80(35, 1.0, Cdchs(i));
end
fprintf('\n-- Discharge C-rate sweep (T=35°C, Cch=1C) --\n');
for i = 1:numel(Cdchs)
    if isnan(life_Cd(i))
        fprintf('Cdch=%.1fC : EFC@80%% > limit (censored)\n', Cdchs(i));
    else
        fprintf('Cdch=%.1fC : EFC@80%% = %.0f\n', Cdchs(i), life_Cd(i));
    end
end
if all(isfinite(life_Cd))
    trendCd = all(diff(life_Cd) <= 0);  % non-increasing vs Cdch
    fprintf('Check (life should decrease with Cdch): %s\n', tern(trendCd,'OK','WARN'));
end


%% 3) Charge C-rate sweep at 35°C (expect weaker effect than Cdch)
Cchs = [2 3 4 5];
life_Cc = nan(size(Cchs));
for i = 1:numel(Cchs)
    life_Cc(i) = efc80(35, Cchs(i), 1.0);
end
fprintf('\n-- Charge C-rate sweep (T=35°C, Cdch=1C) --\n');
for i = 1:numel(Cchs)
    if isnan(life_Cc(i))
        fprintf('Cch=%.1fC : EFC@80%% > limit (censored)\n', Cchs(i));
    else
        fprintf('Cch=%.1fC : EFC@80%% = %.0f\n', Cchs(i), life_Cc(i));
    end
end
if all(isfinite(life_Cc))
    trendCc = all(diff(life_Cc) <= 0);  % non-increasing vs Cch
    fprintf('Check (life should decrease with Cch): %s\n', tern(trendCc,'OK','WARN'));
end

% Optional comparative statement (Cdch effect stronger than Cch)
if all(isfinite([life_Cd life_Cc]))
    weaker = median(abs(diff(life_Cc))) <= 0.7*median(abs(diff(life_Cd)));
    fprintf('\nRelative strength check (Cdch > Cch effect): %s\n', tern(weaker,'OK','WARN'));
end

fprintf('\nDone.\n');



%==================== helpers ====================%

function EFC80 = efc_at_80pct(p, T_degC, Cch, Cdch)
    % Return EFC when SOH first reaches 80%. NaN if not reached.

    % Local params
    pp = p;
    pp.modes.do_calendar = false;
    pp.SOC_min = 0.0;  pp.SOC_max = 1.0;  pp.SOC0 = 0.0;
    pp.T_degC  = T_degC;
    pp.P_ch_max  = 1e12;
    pp.P_dis_max = 1e12;

    % Map C-rates → powers and step durations
    I_ch   =  Cch  * pp.Q_nom;                 % A
    I_dis  =  Cdch * pp.Q_nom;                 % A
    Pch    = (pp.V_nom * I_ch )/1000;          % kW (magnitude)
    Pdis   = (pp.V_nom * I_dis)/1000;          % kW (magnitude)
    h_ch   = 1/max(Cch,  eps);                 % h for full charge (100% DOD)
    h_dis  = 1/max(Cdch, eps);                 % h for full discharge

    % Build many cycles: durations and powers line up 1:1
    Ncyc_max = 12000;                          % generous cap
    dur_steps = repmat([h_ch;  h_dis], Ncyc_max, 1);       % 2*N steps
    Pnet_steps= repmat([-Pch; +Pdis], Ncyc_max, 1);        % 2*N steps

    t = [0; cumsum(dur_steps)];                % time has length = steps+1
    % Split Pnet = Pload - Ppv (use Ppv=0 for clarity), match time length
    Ppv   = zeros(numel(t),1);
    Pload = [Pnet_steps; 0];                   % last value unused by core_loop

    % Run model (time in hours)
    out = core_loop(hours(t), Ppv, Pload, pp);

    % Pick first time SOH ≤ 0.80
    idx = find(out.SOH <= 0.80, 1, 'first');
    if isempty(idx)
        EFC80 = NaN;                           % censored (didn’t reach 80%)
    else
        EFC80 = out.EFC(idx);
    end

end
function s = tern(c, a, b)
% tiny ternary helper: returns a if c is true, else b
if c, s = a; else, s = b; end
end

