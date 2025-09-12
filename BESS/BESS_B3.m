% BESS_B3.m — Energy balance + cycle + calendar aging
clear; clc;

[time, Ppv, Pload] = load_profile_bess();
p = params('B3');                     % cycle ON, calendar ON

out = core_loop(time, Ppv, Pload, p);

plot_bess(out.time, Ppv, Pload, out.SOC, out.SOH);
report(out);
