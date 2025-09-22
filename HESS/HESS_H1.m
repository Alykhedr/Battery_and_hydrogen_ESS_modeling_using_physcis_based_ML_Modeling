
% HESS_M1.m — Baseline HESS (mass balance only, no compressor, no electrochem)
clear; clc;

% 1) load profile
[time, P_el, P_fc] = load_profile_hess();

% 2) params (M1 gates OFF compressor & electrochem)
p = params_hess('M1', time);

% 3) run
out = core_loop_hess(time, P_el, P_fc, p);

% 4) plot & report
plot_hess(out.time, out.P_el, out.P_fc, out.mass_H2, out.m_dot_el, out.m_dot_fc);
report_hess(out);
