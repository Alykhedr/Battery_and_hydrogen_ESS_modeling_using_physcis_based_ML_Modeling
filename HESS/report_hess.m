function report_hess(out)
% REPORT_HESS  Print concise summary of HESS run.

t = out.time;
if isdatetime(t)
    t_h = hours(t - t(1));         % numeric hours since start
elseif isduration(t)
    t_h = hours(t);                 % numeric hours
else
    t_h = t;                        % already hours
end
total_h   = t_h(end) - t_h(1);
span_days = total_h / 24;

prod_kg_total = trapz(t_h, max(out.m_dot_el,0));   % kg (m_dot in kg/h)
cons_kg_total = trapz(t_h, max(out.m_dot_fc,0));   % kg

prod_avg_day = prod_kg_total / span_days;          % kg/day
cons_avg_day = cons_kg_total / span_days;          % kg/day

Ecomp_kWh = trapz(t_h, out.P_comp_kW);             % kWh

fprintf('Final H2 mass = %.2f kg  |  Tank pressure = %.1f bar\n', ...
        out.mass_H2(end), out.P_tank_bar(end));
fprintf('--- SUMMARY ---\n');
fprintf('Span: %.1f days (%.2f years)\n', span_days, span_days/365);
fprintf('Avg production: %.2f kg/day\n', prod_avg_day);
fprintf('Avg consumption: %.2f kg/day\n', cons_avg_day);
fprintf('Compressor energy: %.1f kWh total\n', Ecomp_kWh);
end
