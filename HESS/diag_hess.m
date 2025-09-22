
function diag_hess(out, p)
% Prints quick sanity metrics.

% --- time axis to numeric hours ---
t = out.time;
if isdatetime(t), t_h = hours(t - t(1));
elseif isduration(t), t_h = hours(t);
else, t_h = t - t(1);
end
dt = median(diff(t_h));

% --- aggregates ---
Ecomp_kWh      = trapz(t_h, out.P_comp_kW);
m_prod_kg      = trapz(t_h, max(out.m_dot_el,0));
Ecomp_per_kg   = Ecomp_kWh / max(m_prod_kg, 1e-9);           % kWh/kg
Pel_mean_kW    = mean(max(out.P_el,0));
Pel_stack_kW   = mean(max(out.P_el_remain,0));               % M3 stack power
CF_el_pct      = 100 * Pel_stack_kW / p.P_el_rated;
Pfc_req_mean   = mean(max(out.P_fc,0));
Pfc_impl_mean  = mean(out.m_dot_fc) * p.eta_fc * p.LHV_H2;   % kW
Ptank_minmeanmax= [min(out.P_tank_bar) mean(out.P_tank_bar) max(out.P_tank_bar)];
mass_delta     = out.mass_H2(end) - out.mass_H2(1);
mass_int_check = trapz(t_h, out.m_dot_el - out.m_dot_fc);

% --- print ---
fprintf('dt = %.3f h | span = %.1f h\n', dt, t_h(end)-t_h(1));
fprintf('E_comp = %.1f kWh | E_comp per kg = %.2f kWh/kg\n', Ecomp_kWh, Ecomp_per_kg);
fprintf('PV→EL: mean avail = %.1f kW | mean stack = %.1f kW | CF_EL = %.1f%%\n', ...
        Pel_mean_kW, Pel_stack_kW, CF_el_pct);
fprintf('FC: mean req = %.1f kW | implied from m_dot_fc = %.1f kW\n', ...
        Pfc_req_mean, Pfc_impl_mean);
fprintf('Tank p [bar]: min/mean/max = %.1f / %.1f / %.1f\n', Ptank_minmeanmax);
fprintf('Mass Δ = %.2f kg | ∫(m_el - m_fc)dt = %.2f kg\n', mass_delta, mass_int_check);
end
