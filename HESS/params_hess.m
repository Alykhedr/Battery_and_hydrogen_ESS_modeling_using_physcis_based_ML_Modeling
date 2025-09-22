function p = params_hess(model, time)
% PARAMS_HESS  Parameter set with per-model gating: 'M1' | 'M2' | 'M3'

% normalize tag
if isa(model,'string'); model = char(model); end
model = upper(strtrim(model));
p.meta.model = model;

% ----- Enables by model (match your modelLevel logic) -----
switch model
    case 'M1'
        p.enable.compressor  = false;
        p.enable.electrochem = false;
    case 'M2'
        p.enable.compressor  = true;
        p.enable.electrochem = false;
    case 'M3'
        p.enable.compressor  = true;
        p.enable.electrochem = true;
    otherwise
        error('Unknown model "%s". Use M1 | M2 | M3.', model);
end

% ----- Constants  -----
p.LHV_H2      = 33.33;     % kWh/kg
p.HHV_H2      = 39.35;     % kWh/kg
p.eta_el      = 0.6;       % fixed electrolyzer η (M1/M2 path)
p.eta_fc      = 0.55;      % LHV basis
p.max_mass_H2 = 350;       % kg
p.M_H2        = 0.002016;  % kg/mol
p.tank.T      = 293;       % K
p.tank.V      = 6;         % m^3
p.R           = 8.314;     % J/(mol*K)
p.F           = 96485;     % C/mol
p.P_el_rated  = 1000;      % kW
p.P_fc_rated  = 400;       % kW

% Pre-allocate state arrays
N = numel(time);
p.init.m_dot_el = zeros(1,N);
p.init.m_dot_fc = zeros(1,N);
p.init.mass_H2  = zeros(1,N);
p.init.mass_H2(1) = 100;   % initial tank mass [kg]

% Compressor defaults (only meaningful if enabled)
p.Cp_H2   = 14.3e3 * p.M_H2;  % J/(mol*K)
p.eta_C   = 0.7;
p.gamma   = 1.4;
p.P_inlet = 20  * 1e5;        % Pa
p.P_outlet= 700 * 1e5;        % Pa

% Electrochem defaults (only used if enabled)
p.el.N        = 185;          % cells
p.el.A        = 900;          % cm^2
p.el.V_rev    = 1.23;         % V
p.el.Ohmic    = 0.1;          % ohm*cm^2
p.el.T        = 353;          % K
p.el.P_H2     = 0.54;            % bar
p.el.P_O2     = 0.54;            % bar

% --- Membrane ohmic (Eq. 36–37) ---
p.el.delta_mem_cm   = 175e-4;      % [cm]   δ_mem  (175×10^-4 cm)
p.el.C_Hp_mol_cm3   = 1000/1e6;    % [mol/cm^3]  C_H+ (1000 mol/m^3 → 1e-3 mol/cm^3)
p.el.D_Hp_cm2_s = 2.4e-5;   % [cm^2/s]  D_H+ at Tcell = 60°C  (given: 2.4×10^-9 m^2/s → 2.4×10^-5 cm^2/s)
% NOTE for later: at Tcell = 80°C, D_H+ = 3.0e-5 cm^2/s  (given: 3×10^-9 m^2/s)

% --- Roughness factor inputs (Eq. 34) ---
p.el.phiI_an        = 0.75;        % [-]    φ_I (anode)
p.el.phiI_ca        = 0.75;        % [-]    φ_I (cathode)
p.el.mM_an          = 1.0e-3;      % [g/cm^2]  m_M (anode)
p.el.mM_ca          = 0.3e-3;      % [g/cm^2]  m_M (cathode)
p.el.rhoM_an        = 11.7;         % [g/cm^3]  ρ_M (IrO2)  ← fill from your ref table
p.el.rhoM_ca        =  21.45;       % [g/cm^3]  ρ_M (Pt)
p.el.dM_an          = 2.9e-7;      % [cm]   d_M (anode)  2.9 nm → 2.9e-7 cm
p.el.dM_ca          = 2.7e-7;      % [cm]   d_M (cathode) 2.7 nm → 2.7e-7 cm
p.el.DG0_R          = 237.2e3;
p.z                 = 2;
p.std.P_an_ref_bar = 1.01325; % bar 
p.std.P_cat_ref_bar = 1.01325; % bar
p.el.a_H2O = 1.0; % water activity (liq) ~1

% --- Exchange current Arrhenius (Eq. 35) ---
p.el.Tref           = 298;         % [K]    T_ref (from table)
p.el.i0_ref_an      = 5e-12;         % [A/cm^2] i0_ref (anode) 
p.el.i0_ref_ca      = 1e-2;         % [A/cm^2] i0_ref (cathode) 
p.el.Ea_an          = 76e3;        % [J/mol] Ea (anode)
p.el.Ea_ca          = 4.3e3;       % [J/mol] Ea (cathode)
p.el.L        = 175e-4;            % cm  (membrane thickness)
p.el.alpha_an = 1.2;
p.el.alpha_ca = 0.5;

% Permeability
p.el.epsilon_H2f = 1e-12;     % mol/(cm*s*bar)
p.el.epsilon_O2f = 1e-10;     % mol/(cm*s*bar)
p.el.epsilon_H2p = 0;     % mol/(cm*s*bar)
p.el.beta_H2O    = 1;

% Fixed design current density for compressor sizing/solver bracket
p.el.i        = 5;            % A/cm^2
end
