function p = test_param()
% LFP/graphite 26650 cell — US26650FTC1 (JES 2018)
p.V_nom_V   = 3.2;          % V
p.Q_nom_Ah  = 3.0;          % Ah (BOL)
% current limits (datasheet/test conditions)
p.I_ch_max_A  = 3.0;        % ≈1C charge limit used in the study
p.I_dis_max_A = 20.0;       % continuous discharge limit
p.P_ch_max_kW  = p.V_nom_V*p.I_ch_max_A/1000;
p.P_dis_max_kW = p.V_nom_V*p.I_dis_max_A/1000;

% round-trip energy eff. used only in power reporting
p.eta_en_ch  = 0.99;
p.eta_en_dis = 0.99;

% SOC window to match validation profile
p.SOC_min = 0.054;          % ≈ 5.4%
p.SOC_max = 0.80;           % 80%
p.SOC0    = 0.28;           % 28% start

% environment (setpoint)
p.T_degC  = 25;             % choose test temp; vary per scenario
p.R_gas   = 8.314462618;

% ---- Degradation law (pick one) ----
% A) C/2-only (classic power-law on Ah, A123-style @C/2)
p.deg.mode       = 'c2_only';
p.deg.z          = 0.552;
p.deg.Ea_Jmol    = 31500;
p.deg.B          = 30330;
p.deg.Q_ref_Ah   = 2.0;     % cell-basis in original fits
p.deg.SOH_min    = 0.50;

% B) Unified (rate-aware) — enable by setting p.deg.mode='unified'
p.deg.unified_lnB  = @(C) 1.226*exp(-0.2797*C) + 9.263;
p.deg.unified_Ea   = @(C) 31700 - 370.3*C; % J/mol
end
