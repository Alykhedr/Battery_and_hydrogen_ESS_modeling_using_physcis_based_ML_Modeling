function p = params(model)
% PARAMS  Parameter set with per-model gating: 'B1' | 'B2' | 'B3'

% normalize tag
if isa(model,'string'); model = char(model); end
model = upper(strtrim(model));
p.meta.model = model;

% -------- Cell (Sony LFP US26650FTC1) --------
p.cell.V_nom     = 3.20;   % V
p.cell.Q_nom     = 3.00;   % Ah

% -------- Pack topology --------
p.Ns = 125;              % series  → ~400 V
p.Np = 667;              % parallels → ~2001 Ah

% -------- Derived pack --------
p.V_nom   = p.Ns * p.cell.V_nom;          % ≈ 400 V
p.Q_nom   = p.Np * p.cell.Q_nom;          % ≈ 2001 Ah
p.E_nom   = (p.V_nom * p.Q_nom) / 1000;   % kWh

% Hardware limits (design caps)
p.P_dis_max = min(800, (p.V_nom * (p.Np * 20.0))/1000);  % kW
p.P_ch_max  = min(800, (p.V_nom * (p.Np *  3.0))/1000);  % kW

% Efficiencies
p.eta_ch  = 0.98;
p.eta_dis = 0.98;

% SOC window & initial SOC
p.SOC_min = 0.05; 
p.SOC_max = 0.95; 
p.SOC0    = 0.50;

% Environment (ambient for now)
p.T_degC  = 25;
p.R_gas   = 8.314462618;   % J/mol/K
p.F       = 96485;         % C/mol

% -------- Gating by model --------
switch model
    case 'B1'
        p.modes.do_cycle        = false;
        p.modes.do_calendar     = false;
        p.modes.use_soc_dependence = false;
    case 'B2'
        p.modes.do_cycle        = true;
        p.modes.do_calendar     = false;
        p.modes.use_soc_dependence = false;
    case 'B3'
        p.modes.do_cycle        = true;
        p.modes.do_calendar     = true;
        p.modes.use_soc_dependence = true;
    otherwise
        error('Unknown model "%s". Use B1 | B2 | B3.', model);
end

% -------- Montes cycle-aging (LFP) --------
p.kcyc     = 0.003414;
p.Tref_K   = 298.0;
p.kT       = 5.8755;
p.kDODc    = 0.0046;      % DOD in PERCENT
p.kCch     = 0.1038;
p.kCdch    = 0.296;
p.kmSOC    = 0.0513;
p.a_montes = 0.869;
p.mSOCref  = 0.42;        % 42% → 0.42

% -------- SOH floor --------
p.deg.SOH_min = 0.50;

% -------- Calendar aging (Naumann/Schimpe) --------
p.deg.kcal_ref    = 3.694e-4; % h^(-0.5) at 25°C & ~50% SOC
p.deg.Ea_cal_Jmol = 20592;    % J/mol
p.deg.alpha       = 0.384;
p.deg.k0          = 0.142;
p.deg.Tref_K      = 298.15;   % 25°C
p.deg.Ua_ref_V    = 0.123;    % graphite anode OCV at ~50% SOC (ref)

% OCV mapping (calibrated so Ua(50%)≈0.123 V)
p.deg.OCVmap.xa0   = 0.09;
p.deg.OCVmap.xa100 = 0.78992288;
end
