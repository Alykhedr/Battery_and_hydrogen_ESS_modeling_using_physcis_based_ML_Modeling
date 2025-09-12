%% HESS_2.m  — 2nd-Level HESS Model with compressor-first allocation (no P_req)
clear; clc; close all;

% --- choose level: set to 1 to turn compressor OFF via params() ---
modelLevel = 2;   % <-- set to 1 for "no compressor" mode

% 1) load profile
[time, P_el, P_fc] = load_profiles();

% 2) grab constants & initial state (compressor enabled iff modelLevel>=2)
p = params(time, modelLevel);

% 3) timestep
dt = time(2) - time(1);
N  = numel(time);

% 4) pre-allocate
n_H2       = zeros(1,N);   % mol in tank
P_tank     = zeros(1,N);   % Pa
P_comp     = zeros(1,N);   % W (actual)
m_dot_el   = p.init.m_dot_el;   % kg/h
m_dot_fc   = p.init.m_dot_fc;   % kg/h
P_tank_bar = zeros(1,N);

% --- Debug header ---
fprintf(' hr | Pavail | P_comp_id | P_comp_act | P_remain | m_dot_el | m_dot_fc | P_tank_bar | PR\n');

% 5) initial state
n_H2(1)       = p.init.mass_H2(1) / p.M_H2;
P_tank(1)     = n_H2(1)*p.R*p.tank.T/p.tank.V;
P_tank_bar(1) = P_tank(1)/1e5;
P_inlet_bar = p.P_inlet / 1e5;   

% 6) main simulation loop
for k = 2:N
  %% a) fuel-cell draw (cap by available mass like Model-3)
  if P_fc(k) > 0
    m_fc_req    = P_fc(k)/(p.eta_fc*p.LHV_H2);     % kg/h requested
    m_fc_cap    = max(0, n_H2(k-1)*p.M_H2 / dt);   % kg/h available this step
    m_dot_fc(k) = min(m_fc_req, m_fc_cap);
  else
    m_dot_fc(k) = 0;
  end

  %% b) tank update with LAST step’s flows (Model-3 timing)
  n_H2(k) = n_H2(k-1) + (m_dot_el(k-1) - m_dot_fc(k-1)) * dt / p.M_H2;
  n_H2(k) = min(max(n_H2(k),0), p.max_mass_H2/p.M_H2);

  %% c) tank pressure (ideal gas; same as Model-3 when Z is off)
  P_tank(k)     = n_H2(k)*p.R*p.tank.T/p.tank.V;
  P_tank_bar(k) = P_tank(k)/1e5;

  %% d) available PV
  P_av = max(0, P_el(k));   % kW

  %% e) compressor-first allocation (only if enabled)
  PR = 0; P_comp_id = 0; P_comp_act = 0; P_remain = P_av;   % safe defaults

  if p.enable.compressor && P_av > 0
      % Notional H2 flow at fixed η_el using full PV (legacy estimate)
      mdot_raw = p.eta_el * P_av / p.LHV_H2;           % kg/h
      mdot_mol = mdot_raw / p.M_H2 / 3600;             % mol/s

      PR = max(P_tank_bar(k) / P_inlet_bar, 1.0 + 1e-9);   % dimensionless

      % Ideal compressor demand for that notional flow (for debug/consistency)
      Wc_id     = mdot_mol * p.Cp_H2 * p.tank.T * (PR^((p.gamma-1)/p.gamma)-1) / p.eta_C; % W
      P_comp_id = Wc_id / 1e3;                        % kW

      % Actual compressor takes priority, capped by PV
      P_comp_act = min(P_comp_id, P_av);              % kW
      P_comp(k)  = P_comp_act * 1e3;                  % W

      % Remainder goes to electrolysis
      P_remain = P_av - P_comp_act;                   % kW
  else
      % Compressor disabled or no PV: all PV to electrolysis, compressor = 0
      P_comp(k) = 0;  % W
  end

  %% f) electrolyzer production (fixed η_el model)
  m_dot_el(k) = p.eta_el * P_remain / p.LHV_H2;      % kg/h

  %% g) debug print
  fprintf('%3d | %6.1f | %9.2f | %10.2f | %8.2f | %8.2f | %8.2f | %10.2f | %5.2f\n', ...
          time(k), P_av, P_comp_id, P_comp_act, P_remain, ...
          m_dot_el(k), m_dot_fc(k), P_tank_bar(k), PR);
end

% 7) plot results
mass_H2   = n_H2 * p.M_H2;    % kg
P_comp_kW = P_comp / 1e3;     % kW
plot_hess(time, P_el, P_fc, mass_H2, m_dot_el, m_dot_fc, P_comp_kW, P_tank_bar);

