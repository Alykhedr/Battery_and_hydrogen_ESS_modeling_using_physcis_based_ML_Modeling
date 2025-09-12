% HESS_4.m  — Level 4 HESS Model (Z-enabled tank) + Z-aware pressure cap
% - Z changes pressure readout (real-gas)
% - NEW: tank mass is capped so P_tank <= P_outlet with Z (realistic capacity)
% - Electrolyzer/compressor logic unchanged (priority preserved)
clear; clc; close all;

%% 0) choose level
modelLevel = 4;           

%% 1) load profiles
[time, P_el, P_fc] = load_profiles();

%% 2) constants & flags
p = params(time, modelLevel);   % assumes p.enable.Z = (modelLevel >= 4)

%% 3) timestep
dt = time(2) - time(1);   % hours
N  = numel(time);

%% 4) pre-allocate variables (keep names)
n_H2        = zeros(1,N);
P_tank      = zeros(1,N);    % Pa (real or ideal)
P_tank_bar  = zeros(1,N);    % bar
m_dot_el    = p.init.m_dot_el;
m_dot_fc    = p.init.m_dot_fc;
P_comp      = zeros(1,N);    % W
eta_el      = zeros(1,N);

% keep these preallocated so plots/exports never crash if electrochem is OFF
i_sol    = zeros(1,N);
V_cell   = zeros(1,N);
eta_F    = zeros(1,N);

mass_H2   = zeros(1,N);
m_H2_prev = zeros(1,N);

%% ---- helper: piecewise-constant Z by pressure (bar) ----
Z_of_bar = @(Pbar) ( ...
    (Pbar < 1  ) * 1.0006 + ...
    (Pbar >=1   & Pbar < 10 ) * 1.0059 + ...
    (Pbar >=10  & Pbar < 50 ) * 1.0297 + ...
    (Pbar >=50  & Pbar < 100) * 1.0601 + ...
    (Pbar >=100 & Pbar < 300) * 1.1879 + ...
    (Pbar >=300 & Pbar < 500) * 1.3197 + ...
    (Pbar >=500)              * 1.6454 );

%% ---- NEW: compute Z-aware maximum mass allowed by design pressure ----
if isfield(p,'enable') && isfield(p.enable,'Z') && p.enable.Z
    Pmax      = p.P_outlet;           % Pa (design/target storage pressure)
    Zmax      = Z_of_bar(Pmax/1e5);   % piecewise Z at that pressure
    n_max     = Pmax * p.tank.V / (Zmax * p.R * p.tank.T);  % mol
    m_max_Z   = n_max * p.M_H2;                                % kg
else
    m_max_Z   = inf;  % no Z-cap if Z disabled
end
m_cap = min(p.max_mass_H2, m_max_Z);   % enforce both inventory & pressure limits

%% Debug header
fprintf(' hr | Pavail | P_comp_id | P_comp_act | P_remain | P_el_act | P_bar(real) |  Z | m_cap\n');

%% 5) initial tank values (mass state; pressure via Z/ideal) with Z-cap
mass_H2(1) = min(p.init.mass_H2(1), m_cap);   % <-- cap at start
n_H2(1)    = mass_H2(1)/p.M_H2;               % mol
Pideal     = n_H2(1)*p.R*p.tank.T/p.tank.V;   % Pa
if isfield(p,'enable') && isfield(p.enable,'Z') && p.enable.Z
    Zk          = Z_of_bar(Pideal/1e5);
    P_tank(1)   = Zk * Pideal;
else
    Zk          = 1.0;
    P_tank(1)   = Pideal;
end
P_tank_bar(1) = P_tank(1)/1e5;
fprintf('%3d | %6.1f | %10.2f | %10.2f | %8.2f | %8.2f | %11.2f | %4.2f | %7.2f\n', ...
        time(1), 0, 0, 0, 0, 0, P_tank_bar(1), Zk, m_cap);

%% 6) simulation loop
for k = 2:N
  %% a) fuel-cell consumption (simple mass draw; detailed FC = Level 5)
  if P_fc(k) > 0
    m_fc_req    = P_fc(k)/(p.eta_fc*p.LHV_H2);       % kg/h requested
    m_fc_cap    = max(0, mass_H2(k-1) / dt);         % kg/h available
    m_dot_fc(k) = min(m_fc_req, m_fc_cap);
  else
    m_dot_fc(k) = 0;
  end
  m_H2_prev(k) = mass_H2(k-1);

  %% b) tank update (mass balance) + Z-cap
  mass_H2(k) = mass_H2(k-1) + (m_dot_el(k-1) - m_dot_fc(k-1)) * dt;  % kg
  mass_H2(k) = min(max(mass_H2(k),0), m_cap);                        % <-- Z-cap applied
  n_H2(k)    = mass_H2(k)/p.M_H2;

  % pressure readout (ideal -> Z)
  Pideal = n_H2(k)*p.R*p.tank.T/p.tank.V;   % Pa
  if p.enable.Z
      Zk          = Z_of_bar(Pideal/1e5);
      P_tank(k)   = Zk * Pideal;
  else
      Zk          = 1.0;
      P_tank(k)   = Pideal;
  end
  P_tank_bar(k) = P_tank(k)/1e5;

  %% c) available PV
  P_av = max(0, P_el(k));   % kW

  %% d) compressor-first allocation (unchanged priority)
  P_comp_id = 0; P_comp_act = 0; P_remain = P_av;  % defaults

  if p.enable.compressor
    if p.enable.electrochem
      % detailed electrolyzer informs max molar flow to compress
      V_design   = cell_voltage(p.el.i, p);
      area_cm2   = p.el.A * 1e4;
      P_req      = p.el.N * p.el.i * area_cm2 * V_design / 1e3;     % kW
      frac_max   = min(1, P_av / max(P_req, eps));
      i_max      = p.el.i * frac_max;
      V_max      = cell_voltage(i_max, p);
      I_tot_max  = p.el.N * i_max * area_cm2;
      eta_fmax   = faraday_eff(i_max, p);
      mol_s_max  = I_tot_max * eta_fmax / (2*p.F);                  % mol/s

      % fixed setpoint PR path (keep behavior consistent with your Model-3)
      PR         = max(p.P_outlet / p.P_inlet, 1.0);
      Wc_id      = mol_s_max * p.Cp_H2 * p.tank.T * (PR^((p.gamma-1)/p.gamma)-1) / p.eta_C; % W
      P_comp_id  = Wc_id / 1e3;                                     % kW
      P_comp_act = min(P_comp_id, P_av);                            % kW
      P_comp(k)  = P_comp_act * 1e3;                                % W
      P_remain   = P_av - P_comp_act;                               % kW
    else
      % fixed-η electrolyzer path (Model-2 style)
      mdot_raw   = p.eta_el * P_av / p.LHV_H2;                      % kg/h
      mdot_mol   = mdot_raw / p.M_H2 / 3600;                        % mol/s
      PR         = max(P_tank(k) / p.P_inlet, 1.0);
      Wc_id      = mdot_mol * p.Cp_H2 * p.tank.T * (PR^((p.gamma-1)/p.gamma)-1) / p.eta_C; % W
      P_comp_id  = Wc_id / 1e3;                                     % kW
      P_comp_act = min(P_comp_id, P_av);                            % kW
      P_comp(k)  = P_comp_act * 1e3;                                % W
      P_remain   = P_av - P_comp_act;                               % kW
    end
  else
    P_comp(k) = 0;  % compressor off
  end

  %% e) electrolyzer power & production
  if p.enable.electrochem
    if P_remain > 0
      V_design = cell_voltage(p.el.i, p);
      area_cm2 = p.el.A * 1e4;
      P_req    = p.el.N * p.el.i * area_cm2 * V_design / 1e3;  % kW

      if P_remain >= P_req
        i_act = p.el.i;
      else
        fun  = @(i) p.el.N * i * area_cm2 .* cell_voltage(i,p) / 1e3 - P_remain;
        i_act = fzero(fun, [0, p.el.i]);
      end
      V_act        = cell_voltage(i_act, p);
      P_el_act     = p.el.N * i_act * area_cm2 * V_act / 1e3;          % kW
      eta_F(k)     = faraday_eff(i_act, p);
      mol_s        = p.el.N * i_act * area_cm2 * eta_F(k) / (2*p.F);   % mol/s
      m_dot_el(k)  = mol_s * p.M_H2 * 3600;                            % kg/h
      eta_el(k)    = (P_el_act>0) * (m_dot_el(k) * p.LHV_H2 / max(P_el_act,eps));
      i_sol(k)     = i_act; 
      V_cell(k)    = V_act;
    else
      P_el_act     = 0; 
      m_dot_el(k)  = 0;
    end
  else
    % fixed-η electrolyzer (Model-2 behavior)
    P_el_act     = P_remain;                                          % kW
    m_dot_el(k)  = p.eta_el * P_el_act / p.LHV_H2;                    % kg/h
    eta_el(k)    = (P_el_act>0) * (m_dot_el(k) * p.LHV_H2 / max(P_el_act,eps));
  end

  %% f) debug print
  fprintf('%3d | %6.1f | %10.2f | %10.2f | %8.2f | %8.2f | %11.2f | %4.2f | %7.2f\n', ...
          time(k), P_av, P_comp_id, P_comp_act, P_remain, P_el_act, P_tank_bar(k), Zk, m_cap);
end

%% 7) plots
figure
plot(i_sol, eta_el, '.-','LineWidth',1.5)
xlabel('Current density i (A/cm^2)'), ylabel('\eta_{el}')
title('Electrolyzer efficiency vs. current density'), grid on

figure
subplot(3,1,1); plot(time, V_cell, 'LineWidth',1.5); ylabel('V_{cell} (V)'); title('Electrolyzer Cell Voltage'); grid on
subplot(3,1,2); plot(time, i_sol,  'LineWidth',1.5); ylabel('i (A/cm^2)');   title('Actual Current Density');   grid on
subplot(3,1,3); plot(time, eta_F,  'LineWidth',1.5); ylabel('\eta_F'); xlabel('Time (h)'); title('Faraday Efficiency over Time'); grid on

%% 8) full system plot
P_comp_kW = P_comp/1e3;
plot_hess(time, P_el, P_fc, mass_H2, m_dot_el, m_dot_fc, P_comp_kW, P_tank_bar);

%% 9) export
tbl = table(time(:), P_el(:), P_fc(:), m_H2_prev(:), m_dot_el(:), m_dot_fc(:), eta_el(:), ...
            i_sol(:), V_cell(:), eta_F(:), P_comp_kW(:), mass_H2(:), P_tank_bar(:), ...
            'VariableNames', {'Hour','P_el','P_fc','m_H2_prev','m_dot_el','m_dot_fc','eta_el', ...
                              'i_sol','V_cell','eta_F','P_comp','m_H2','P_tank_bar'});
writetable(tbl, 'model4_train_data.csv');
fprintf('Exported %d rows\n', height(tbl));
