
% HESS_M3.m — Electrochemical + compressor-first allocation
clear; clc;

% 1) load profile
[time, P_el, P_fc] = load_profile_hess();

% 2) params (M3: compressor ON, electrochem ON)
p = params_hess('M3', time);

% 3) run
out = core_loop_hess(time, P_el, P_fc, p);

% 4) plots (electrochem panels shown automatically if enabled)
plot_hess(out.time, out.P_el, out.P_fc, out.mass_H2, out.m_dot_el, out.m_dot_fc, ...
          out.P_comp_kW, out.P_tank_bar, out.i_sol, out.V_cell, out.eta_F);
report_hess(out);

