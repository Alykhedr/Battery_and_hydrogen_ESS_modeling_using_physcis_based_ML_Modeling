% params.m
function p = params(time, modelLevel)

% Determine which model to enable
p.modelLevel         = modelLevel;
p.enable.compressor  = (modelLevel >= 2);
p.enable.electrochem = (modelLevel >= 3);
p.enable.Z           = (modelLevel >= 4);  
p.enable.fuelcell    = (modelLevel >= 5);

% Constants
p.LHV_H2      = 33.33;    % kWh/kg
p.HHV_H2      = 39.35;    % kWh/kg
p.eta_el      = 0.6;      % 60% LHV
p.eta_fc      = 0.55;     % 49% LHV
p.max_mass_H2 = 350;      % kg
p.M_H2        = 0.002016; % kg/mol
p.tank.T      = 293;      % K
p.tank.V      = 6;        % m^3
p.R           = 8.314;    % J/(mol*K)
p.F           = 96485;    % C/mol

% Rated sizes
p.P_el_rated = 1000;  % kW
p.P_fc_rated = 400;   % kW

% Pre-allocate state arrays
N = numel(time);
p.init.m_dot_el = zeros(1,N);
p.init.m_dot_fc = zeros(1,N);
p.init.mass_H2  = zeros(1,N);

% Initial tank level
p.init.mass_H2(1) = 100;

% Compressor defaults
if p.enable.compressor
  p.Cp_H2    = 14.3e3 * p.M_H2;  % J/(mol*K)
  p.eta_C    = 0.7;
  p.gamma    = 1.4;
  p.P_inlet  = 20   * 1e5;       % Pa
  p.P_outlet = 350  * 1e5;       % Pa

end

% Electrochemistry defaults
if p.enable.electrochem
  p.el.N        = 185;           % cells for 1MW @3A/cm2 & 2V
  p.el.A        = 900 * 1e-4;    % m^2 (30x30 cm)
  p.el.V_rev       = 1.23;          % V
  p.el.Ohmic       = 0.1;          % ohm*cm^2
  p.el.T        = 370;           % K
  p.el.P_H2        = 5;             % bar
  p.el.P_O2        = 5;             % bar

  % Kinetic & transport
  p.el.i_p         = 1e-3;          % A/cm^2
  p.el.i_n         = 1e-3;          % A/cm^2
  p.el.sigma_p     = 0.5;
  p.el.sigma_n     = 0.5;
  p.el.L           = 0.05;          % cm

  % Permeability
  p.el.epsilon_H2f = 1e-12;         % mol/(cm*s*bar)
  p.el.epsilon_O2f = 1e-10;         % mol/(cm*s*bar)
  p.el.epsilon_H2p = 1e-12;         % mol/(cm*s*bar)
  p.el.beta_H2O    = 1;

  % Fixed design current density
  p.el.i  = 3;                      % A/cm^2
end

if p.enable.fuelcell
      p.FC.R0    = 0.1537;      % Eq (3) pre‑exp factor [Ω]
      p.FC.N = 100;   % ← number of cells in the stack
    p.FC.Ea_R  = 1800;        % Eq (3) activation energy [J/mol]
    p.FC.A0    = 0.1591;      % Eq (4) pre‑exp Tafel slope [V]
    p.FC.Ea_A  = 5344;        % Eq (4) activation energy [J/mol]
    p.FC.Iex   = 1e-6;        % Exchange current [A]
    p.FC.P_H2  = 1.35;        % Eq (2) H2 partial pressure [atm]
    p.FC.P_O2  = 1.00;        % Eq (2) O2 partial pressure [atm]
    % (thermal dynamics)
    p.FC.T1    = 353;         % Eq (5) initial stack T [K]
    p.FC.T2    = 300;         % Eq (5) final/ambient T [K]
    p.FC.Ht    = 15.07;       % Eq (5) heat‑transfer coeff [W/K]
    p.FC.mcp   = 4304;        % Eq (5) heat capacitance [J/K]
end


end