% HESS_4_clean.m
% Level 4 HESS Model with fixed FC current (no root‐finder)
clear; clc; close all;

%% 1) load profiles
[time,P_el,P_fc] = load_profiles();

%% 2) parameters
p = params(time,4);

%% 3) time discretization
dt = time(2) - time(1);
N  = numel(time);

%% 4) pre‑allocate all time‐series
n_H2       = zeros(1,N);
mass_H2    = zeros(1,N);
P_tank     = zeros(1,N);
P_tank_bar = zeros(1,N);

m_dot_el   = zeros(1,N);
m_dot_fc   = zeros(1,N);
I_FC       = zeros(1,N);

P_comp     = zeros(1,N);
Pcid       = zeros(1,N);
Pc_act     = zeros(1,N);
P_rem      = zeros(1,N);
P_el_act   = zeros(1,N);

i_el       = zeros(1,N);
eta_el     = zeros(1,N);
eta_F      = zeros(1,N);

%% cell areas [cm^2]
area_el_cm2 = p.el.A * 1e4;
if isfield(p.FC,'cell_area')
    area_fc_cm2 = p.FC.cell_area * 1e4;
elseif isfield(p.FC,'area')
    area_fc_cm2 = p.FC.area * 1e4;
else
    warning('No FC cell‐area in params, using electrolyzer area.');
    area_fc_cm2 = area_el_cm2;
end

%% 5) nominal FC stack voltage
V_cell_nom  = 0.8;              % per‐cell voltage [V]
V_stack_nom = V_cell_nom * p.FC.N;

%% 6) initialize tank
n_H2(1)       = p.init.mass_H2(1)/p.M_H2;
mass_H2(1)    = p.init.mass_H2(1);
P_tank(1)     = n_H2(1)*p.R*p.tank.T/p.tank.V;
P_tank_bar(1) = P_tank(1)/1e5;

fprintf(' hr | Pavail | Pcid   | Pc_act | Premain | Pel_act | P_tank_bar\n');

%% 7) main loop
for k = 2:N
  %% a) FC: fixed V, so I = P/V_stack_nom
  if P_fc(k)>0 && mass_H2(k-1)>0
    I_FC(k)     = P_fc(k)*1e3 / V_stack_nom;           % A
    m_dot_fc(k) = I_FC(k)/(2*p.F)*p.M_H2*3600;         % kg/h H2
    eta_F(k)    = (V_stack_nom*I_FC(k)/1e3) / P_fc(k);
  else
    I_FC(k)     = 0;
    m_dot_fc(k) = 0;
    eta_F(k)    = 0;
  end

  %% b) tank update
  n_H2(k)       = n_H2(k-1) + (m_dot_el(k-1)-m_dot_fc(k-1))*dt/p.M_H2;
  n_H2(k)       = min(max(n_H2(k),0), p.max_mass_H2/p.M_H2);
  mass_H2(k)    = n_H2(k)*p.M_H2;
  P_tank(k)     = n_H2(k)*p.R*p.tank.T/p.tank.V;
  P_tank_bar(k) = P_tank(k)/1e5;

  %% c) compressor‑first
  P_avail = max(0, P_el(k));
  V_design= cell_voltage(p.el.i,p);
  P_req   = p.el.N*p.el.i*area_el_cm2*V_design/1e3;

  if p.enable.compressor
    frac       = min(1, P_avail/P_req);
    i_max      = p.el.i * frac;
    V_max      = cell_voltage(i_max,p);
    I_tot      = p.el.N * i_max * area_el_cm2;
    eta_f      = faraday_eff(i_max,p);
    mol_s      = I_tot*eta_f/(2*p.F);
    PR         = p.P_outlet / max(P_tank(k), p.P_inlet);
    Wc         = mol_s * p.Cp_H2 * p.tank.T * (PR^((p.gamma-1)/p.gamma)-1) / p.eta_C;
    Pcid(k)    = Wc/1e3;
    Pc_act(k)  = min(Pcid(k), P_avail);
    P_comp(k)  = Pc_act(k)*1e3;
    P_rem(k)   = P_avail - Pc_act(k);
  else
    Pcid(k)   = 0;
    Pc_act(k) = 0;
    P_comp(k) = 0;
    P_rem(k)  = P_avail;
  end

  %% d) electrolyzer on the remainder
  if P_rem(k)>0
    if P_rem(k) >= P_req
      i_act = p.el.i;
    else
      fun_el = @(i) p.el.N*i*area_el_cm2 .* cell_voltage(i,p)/1e3 - P_rem(k);
      i_act  = fzero(fun_el, [0 p.el.i]);
    end
    V_act       = cell_voltage(i_act,p);
    P_el_act(k) = p.el.N*i_act*area_el_cm2 * V_act/1e3;
    m_dot_el(k) = faraday_eff(i_act,p)*p.el.N*i_act*area_el_cm2/(2*p.F)*p.M_H2*3600;
    eta_el(k)   = m_dot_el(k)*p.LHV_H2 / P_el_act(k);
    i_el(k)     = i_act;
  else
    P_el_act(k) = 0;
    m_dot_el(k) = 0;
    i_el(k)     = 0;
    eta_el(k)   = 0;
  end

  %% e) print debug line
  fprintf('%3d | %6.1f | %6.2f | %6.2f | %8.2f | %7.2f | %10.2f\n', ...
          time(k), P_avail, Pcid(k), Pc_act(k), P_rem(k), P_el_act(k), P_tank_bar(k));
end

%% 8) Plot FC & electrolyzer currents / densities
figure;
subplot(2,1,1)
plot(time, I_FC, '-','LineWidth',1.5)
ylabel('I_{FC} (A)'); title('Fuel‐Cell Stack Current'); grid on

subplot(2,1,2)
j_FC = I_FC / area_fc_cm2;  % A/cm^2 per cell
plot(time, j_FC, '-','LineWidth',1.5)
ylabel('j_{FC} (A/cm^2)'); xlabel('Time (h)')
title('Fuel Cell Current Density'); grid on

figure;
plot(time, i_el, '-','LineWidth',1.5)
xlabel('Time (h)'); ylabel('i_{EL} (A/cm^2)')
title('Electrolyzer Current Density'); grid on
