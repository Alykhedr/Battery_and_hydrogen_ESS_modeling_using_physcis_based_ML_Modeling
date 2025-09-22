
% HESS_M2.m — Compressor-first allocation (compressor ON, electrochem OFF)
clear; clc;

% 1) load profile
[time, P_el, P_fc] = load_profile_hess();

% 2) params (M2: compressor ON, electrochem OFF)
p = params_hess('M2', time);

% 3) run
out = core_loop_hess(time, P_el, P_fc, p);

% 4) plot & report
plot_hess(out.time, out.P_el, out.P_fc, out.mass_H2, out.m_dot_el, out.m_dot_fc, ...
          out.P_comp_kW, out.P_tank_bar);
report_hess(out);
