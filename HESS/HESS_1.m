% HESS Model 1  
% Author: Aly Khedr

clear; clc; close all;

% 1) load profile
[time, P_el, P_fc] = load_profiles();

% 2) grab constants + initial state
p        = params(time,1);
p.enable.compressor = false; 
m_dot_el = p.init.m_dot_el;
m_dot_fc = p.init.m_dot_fc;
mass_H2  = p.init.mass_H2;

% 3) compute dt
dt = time(2) - time(1);

% 4) run the HESS loop
for k = 2:numel(time)
  if P_el(k)>0
    m_dot_el(k) = p.eta_el * P_el(k) / p.LHV_H2;
  end
  if P_fc(k)>0
    m_dot_fc(k) = P_fc(k) / (p.eta_fc * p.LHV_H2);
  end

 mass_H2(k) = mass_H2(k-1) + (m_dot_el(k-1) - m_dot_fc(k-1))*dt;
  mass_H2(k) = min(max(mass_H2(k),0), p.max_mass_H2);
end

% 5) plot
plot_hess(time, P_el, P_fc, mass_H2,m_dot_el,m_dot_fc);
