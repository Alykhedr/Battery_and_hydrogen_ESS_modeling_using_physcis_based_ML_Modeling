% run_models.m — one switch to run M1 (mass), M2 (+compressor), or M3 (+electrochem)
clear; clc;

%% ---------------- Choose model here ----------------
% Options: 'M1' (mass balance only), 'M2' (compressor-first), 'M3' (electrochem + compressor)
model = 'M3';

%% --------------- Load profile ----------------------
[time, P_el, P_fc] = yearly_power_profile();
time = time(:); P_el = P_el(:); P_fc = P_fc(:);

%% --------------- Params --------------
p = params_hess(model, time);

%% --------------- Run the unified core --------------
out = core_loop_hess(time, P_el, P_fc, p);

%% --------------- Plots -----------------------------
if p.enable.electrochem
    plot_hess(out.time, out.P_el, out.P_fc, out.mass_H2, out.m_dot_el, out.m_dot_fc, ...
              out.P_comp_kW, out.P_tank_bar);
elseif p.enable.compressor
    plot_hess(out.time, out.P_el, out.P_fc, out.mass_H2, out.m_dot_el, out.m_dot_fc, ...
              out.P_comp_kW, out.P_tank_bar);
else
    plot_hess(out.time, out.P_el, out.P_fc, out.mass_H2, out.m_dot_el, out.m_dot_fc);
end
diag_hess(out,p);

fprintf('PR>1 hours = %.1f%%\n', 100*mean(out.P_tank_bar > p.P_inlet/1e5));

t = out.time; if isdatetime(t), t = hours(t - t(1)); end
mask = out.P_tank_bar > p.P_inlet/1e5;
Ecomp_PR1 = trapz(t(mask), out.P_comp_kW(mask));
mprod_PR1 = trapz(t(mask), max(out.m_dot_el(mask),0));
fprintf('E_comp(PR>1) = %.2f kWh/kg\n', Ecomp_PR1/max(mprod_PR1,1e-9));

%% --------------- Report ----------------------------
report_hess(out);

%% --------------- Export CSV for ML -----------------
if strcmpi(model,'M3')
    tbl = table( ...
        out.time(:), ...
        out.P_el(:), ...
        out.P_fc(:), ...
        out.m_dot_el(:), ...
        out.m_dot_fc(:), ...
        out.eta_el(:), ...
        out.i_model(:), ...
        out.V(:), ...
        out.eta_F(:), ...
        out.P_comp_kW(:), ...
        out.mass_H2(:), ...
        out.P_tank_bar(:), ...
        'VariableNames', {'Hour','P_el','P_fc','m_dot_el','m_dot_fc', ...
                          'eta_el','i_model','V','eta_F','P_comp','m_H2','P_tank_bar'} ...
    );
    writetable(tbl, 'model3_train_data.csv');
    fprintf('Exported %d rows to model3_train_data.csv\n', height(tbl));
end
