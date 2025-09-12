% BESS_B1.m — Baseline energy balance (SOC only)
clear; clc;

% --- Load profile
[time, Ppv, Pload] = load_profile_bess();

% --- Params (B1 gates OFF cycle/calendar)
p = params('B1');

% --- Run
out = core_loop(time, Ppv, Pload, p);

% --- Plot & report
plot_bess(out.time, Ppv, Pload, out.SOC, ones(numel(out.time),1)); % no SOH in B1, pass ones
report(out);
