% BESS_B2.m — Energy balance + cycle aging (Montes)
clear; clc;

[time, Ppv, Pload] = load_profile_bess();
p = params('B2');                     % cycle ON, calendar OFF

out = core_loop(time, Ppv, Pload, p);

plot_bess(out.time, Ppv, Pload, out.SOC, out.SOH);
report(out);
