% HESS_3.m  
% Level 3 HESS Model (Electrochemical + Compressor-first allocation)
clear; clc; close all;

%% 0) choose level
modelLevel = 3;            % <-- set to 2 or 1 to downgrade features

%% 1) load profiles
[time, P_el, P_fc] = load_profiles();

%% 2) constants & flags (drive enables from modelLevel)
p = params(time, modelLevel);  % sets p.enable.compressor / p.enable.electrochem

%% 3) timestep
dt = time(2) - time(1);
N  = numel(time);

%% 4) pre-allocate variables (names unchanged)
n_H2        = zeros(1,N);
P_tank      = zeros(1,N);
P_comp      = zeros(1,N);
eta_el      = zeros(1,N);
P_tank_bar  = zeros(1,N);

% Always allocate these so export/plots don't crash; fill only when used
i_sol       = zeros(1,N);
V_cell      = zeros(1,N);
eta_F       = zeros(1,N);

% time-invariant area only needed if electrochem is on
if p.enable.electrochem
  area_cm2 = p.el.A * 1e4;
end

m_dot_el    = p.init.m_dot_el;
m_dot_fc    = p.init.m_dot_fc;
mass_H2     = zeros(1,N);
m_H2_prev   = zeros(1,N);

%% Debug header
fprintf(' hr | Pavail | P_comp_id | P_comp_act | P_remain | P_el_act | P_tank_bar\n');

%% 5) initial tank values (ideal gas here)
n_H2(1)        = p.init.mass_H2(1)/p.M_H2;
mass_H2(1)     = p.init.mass_H2(1);
P_tank(1)      = n_H2(1)*p.R*p.tank.T/p.tank.V;
P_tank_bar(1)  = P_tank(1)/1e5;
P_inlet_bar = p.P_inlet / 1e5;  


%% 6) simulation loop
for k = 2:N
  % a) fuel-cell consumption (cap by available tank mass like Model-2/3)
  if P_fc(k) > 0
    m_fc_req    = P_fc(k)/(p.eta_fc*p.LHV_H2);   % kg/h
    m_fc_cap    = max(0, mass_H2(k-1)/dt);       % kg/h available this step
    m_dot_fc(k) = min(m_fc_req, m_fc_cap);
  else
    m_dot_fc(k) = 0;
  end
  m_H2_prev(k) = mass_H2(k-1);

  % b) tank update (mass balance using last step’s actual flows)
  n_H2(k)    = n_H2(k-1) + (m_dot_el(k-1) - m_dot_fc(k-1))*dt/p.M_H2;
  n_H2(k)    = min(max(n_H2(k),0), p.max_mass_H2/p.M_H2);
  mass_H2(k) = n_H2(k)*p.M_H2;

  % pressure (ideal gas; your Z work can be re-inserted later if desired)
  P_tank(k)     = n_H2(k)*p.R*p.tank.T/p.tank.V;
  P_tank_bar(k) = P_tank(k)/1e5;

  % c) available PV
  P_available = max(0, P_el(k));  % kW

  % d) compressor-first allocation
  P_comp_id  = 0; P_comp_act = 0; P_remain = P_available; PR = 0;

  if p.enable.compressor && P_available > 0
    % --- Compressor logic depends on whether electrochemistry is active ---
    if p.enable.electrochem
      % (Model-3 path) size compressor for the "i_max" electrolyzer guess
      V_design   = cell_voltage(p.el.i, p);
      P_req      = p.el.N * p.el.i * area_cm2 * V_design / 1e3;      % kW
      frac_max   = min(1, P_available / P_req);
      i_max      = p.el.i * frac_max;
      V_max      = cell_voltage(i_max, p);
      I_tot_max  = p.el.N * i_max * area_cm2;
      eta_fmax   = faraday_eff(i_max, p);
      mol_s_max  = I_tot_max * eta_fmax / (2*p.F);                   % mol/s

      PR = max(P_tank_bar(k) / P_inlet_bar, 1.0 + 1e-9);   % dimensionless
      W_comp     = mol_s_max * p.Cp_H2 * p.tank.T ...
                   * (PR^((p.gamma-1)/p.gamma) - 1) / p.eta_C;       % W
      P_comp_id  = W_comp / 1e3;                                     % kW
      P_comp_act = min(P_comp_id, P_available);                       % kW
      P_comp(k)  = P_comp_act * 1e3;                                  % W
      P_remain   = P_available - P_comp_act;                          % kW

    else
      % (Model-2 behavior) fixed-η electrolyzer estimate for compressor sizing
      mdot_raw   = p.eta_el * P_available / p.LHV_H2;                 % kg/h
      mdot_mol   = mdot_raw / p.M_H2 / 3600;                          % mol/s
      PR         = max(P_tank(k) / p.P_inlet, 1.0);
      W_comp     = mdot_mol * p.Cp_H2 * p.tank.T ...
                   * (PR^((p.gamma-1)/p.gamma) - 1) / p.eta_C;        % W
      P_comp_id  = W_comp / 1e3;                                      % kW
      P_comp_act = min(P_comp_id, P_available);                        % kW
      P_comp(k)  = P_comp_act * 1e3;                                  % W
      P_remain   = P_available - P_comp_act;                          % kW
    end
  else
    % compressor disabled or no PV
    P_comp(k) = 0;   % W
  end

  %% e) electrolyzer block
  if p.enable.electrochem
    % (Model-3) solve exact current for P_remain
    V_design  = cell_voltage(p.el.i, p);
    P_req     = p.el.N * p.el.i * area_cm2 * V_design / 1e3;   % kW

    if P_remain > 0
      if P_remain >= P_req
        i_act = p.el.i;
      else
        fun  = @(i) p.el.N * i * area_cm2 .* cell_voltage(i,p) / 1e3 - P_remain;
        i_act = fzero(fun, [0, p.el.i]);
      end
      V_act      = cell_voltage(i_act, p);
      P_el_act   = p.el.N * i_act * area_cm2 * V_act / 1e3;           % kW
      eta_F(k)   = faraday_eff(i_act, p);
      mol_s      = p.el.N * i_act * area_cm2 * eta_F(k) / (2*p.F);    % mol/s
      m_dot_el(k)= mol_s * p.M_H2 * 3600;                              % kg/h
      eta_el(k)  = m_dot_el(k) * p.LHV_H2 / max(P_el_act, eps);
      i_sol(k)   = i_act;
      V_cell(k)  = V_act;
    else
      P_el_act   = 0;
      m_dot_el(k)= 0;
    end

  else
    % (Model-1/2) fixed-η electrolyzer using the remainder
    P_el_act    = P_remain;                                           % kW
    m_dot_el(k) = p.eta_el * P_el_act / p.LHV_H2;                     % kg/h
    eta_el(k)   = (P_el_act>0) * p.eta_el;                             % for plotting/export
    % i_sol, V_cell, eta_F remain zero
  end

  % f) debug print
  fprintf('%3d | %6.1f | %10.2f | %10.2f | %8.2f | %8.2f | %10.2f\n', ...
          time(k), P_available, P_comp_id, P_comp_act, P_remain, P_el_act, P_tank_bar(k));
end

%% 7) plot overall electrolyzer efficiency vs current density (only if electrochem on)
if p.enable.electrochem
  figure
  plot(i_sol, eta_el, '.-','LineWidth',1.5)
  xlabel('Current density i (A/cm^2)')
  ylabel('\eta_{el}')
  title('Electrolyzer efficiency vs. current density')
  grid on
end

%% 8) plot cell voltage, current & efficiency over time (only if electrochem on)
if p.enable.electrochem
  figure
  subplot(3,1,1)
  plot(time, V_cell, 'LineWidth',1.5)
  ylabel('V_{cell} (V)')
  title('Electrolyzer Cell Voltage')
  grid on

  subplot(3,1,2)
  plot(time, i_sol, 'LineWidth',1.5)
  ylabel('i (A/cm^2)')
  title('Actual Current Density')
  grid on

  subplot(3,1,3)
  plot(time, eta_F, 'LineWidth',1.5)
  ylabel('\eta_F')
  xlabel('Time (h)')
  title('Faraday Efficiency over Time')
  grid on
end

%% 9) full system plot (always)
P_comp_kW = P_comp/1e3;
plot_hess(time, P_el, P_fc, mass_H2, m_dot_el, m_dot_fc, P_comp_kW, P_tank_bar);

%% 10) export
tbl = table(time(:), P_el(:), P_fc(:), m_H2_prev(:), m_dot_el(:), m_dot_fc(:), eta_el(:), ...
            i_sol(:), V_cell(:), eta_F(:), P_comp_kW(:), mass_H2(:), P_tank_bar(:), ...
            'VariableNames', {'Hour','P_el','P_fc','m_H2_prev','m_dot_el','m_dot_fc','eta_el', ...
                              'i_sol','V_cell','eta_F','P_comp','m_H2','P_tank_bar'});
writetable(tbl, 'model3_train_data.csv');
fprintf('Exported %d rows\n', height(tbl));