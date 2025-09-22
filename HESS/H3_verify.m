%% HESS_system_verify.m  — system-level checks for M1/M2/M3 + eta_F
clear; clc; close all;

% ---- Common profile ----
[time, P_el, P_fc] = load_profile_hess();
dt = time(2)-time(1);

% ---- Params ----
pm1 = params_hess('M1', time);
pm2 = params_hess('M2', time);
pm3 = params_hess('M3', time);

% Paper-like stack operating point for M3
pm3.el.T    = 333;                 % 80 C
pm3.el.P_H2 = 0.54; pm3.el.P_O2 = 0.54; pm3.el.a_H2O = 1.0;
% Optional (your calibrated values)
% pm3.el.D_Hp_cm2_s = 3.5e-5;
% pm3.el.i0_ref_ca  = 1e-2;

% ---- Run models ----
out1 = core_loop_hess(time, P_el, P_fc, pm1);
out2 = core_loop_hess(time, P_el, P_fc, pm2);
out3 = core_loop_hess(time, P_el, P_fc, pm3);

% ---- MASS/ENERGY CHECKS ----
fprintf('\n== MASS/ENERGY CHECKS ==\n');
sumstats = @(o) struct( ...
    'dm',   o.mass_H2(end)-o.mass_H2(1), ...
    'prod', trapz(o.time, max(o.m_dot_el,0)) * (24 / (o.time(end)-o.time(1))), ... % kg/day (avg over span)
    'cons', trapz(o.time, max(o.m_dot_fc,0)) * (24 / (o.time(end)-o.time(1))), ... % kg/day (avg over span)
    'Ecmp', trapz(o.time, max(o.P_comp_kW,0)) );                                   % kWh total

S1 = sumstats(out1); S2 = sumstats(out2); S3 = sumstats(out3);
fprintf('M1  Δm=%6.2f kg | prod=%5.2f kg/d | cons=%5.2f kg/d | Ecomp=%6.1f kWh\n', S1.dm,S1.prod,S1.cons,S1.Ecmp);
fprintf('M2  Δm=%6.2f kg | prod=%5.2f kg/d | cons=%5.2f kg/d | Ecomp=%6.1f kWh\n', S2.dm,S2.prod,S2.cons,S2.Ecmp);
fprintf('M3  Δm=%6.2f kg | prod=%5.2f kg/d | cons=%5.2f kg/d | Ecomp=%6.1f kWh\n', S3.dm,S3.prod,S3.cons,S3.Ecmp);

%% Electrolyzer SEC (kWh/kg) — mass-weighted mean only
outs = struct('M2', out2, 'M3', out3);
tags = fieldnames(outs);
for t = 1:numel(tags)
    o = outs.(tags{t});
    on = o.m_dot_el > 1e-9;
    if any(on)
        SEC_mean = sum(o.P_el_remain(on)) / sum(o.m_dot_el(on));   % kWh/kg
        fprintf('\n%s EL-SEC (mean, mass-weighted): %.2f kWh/kg\n', tags{t}, SEC_mean);
    else
        fprintf('\n%s EL-SEC: (no H2 production)\n', tags{t});
    end
end

%% Compressor SEC — mean only (+ DOE in-band share)
DOE_min = 1.7; DOE_max = 6.4; eps_mdot = 1e-9;
outs = struct('M2', out2, 'M3', out3);
tags = fieldnames(outs);
for t = 1:numel(tags)
    o = outs.(tags{t});
    on = (o.P_comp_kW > 1e-9) & (o.m_dot_el > eps_mdot);
    if any(on)
        SECc = o.P_comp_kW(on) ./ o.m_dot_el(on);                   % kWh/kg
        inband = mean(SECc >= DOE_min & SECc <= DOE_max) * 100;
        fprintf('%s COMP-SEC (mean): %.2f kWh/kg  | DOE in-band: %.1f%%\n', ...
                tags{t}, mean(SECc), inband);
    else
        fprintf('%s COMP-SEC: (compressor never on)\n', tags{t});
    end
end

%% Fuel-cell yield (kWh/kg) — mean only (+ limited-steps flag)
eps_m = 1e-9;
runs = { {'M1', out1, pm1}, {'M2', out2, pm2}, {'M3', out3, pm3} };
for r = 1:numel(runs)
    tag = runs{r}{1}; o = runs{r}{2}; p = runs{r}{3};
    on  = (o.P_fc > 0) & (o.m_dot_fc > eps_m);
    if any(on)
        % mass-limited indicator (requested vs delivered)
        m_req = o.P_fc(on) ./ (p.eta_fc * p.LHV_H2);
        limited_pct = 100 * mean(o.m_dot_fc(on) < (m_req - 1e-9));

        SEC_fc_mean = sum(o.P_fc(on)) / sum(o.m_dot_fc(on));       % kWh/kg
        fprintf('%s FC-yield (mean): %.2f kWh/kg  | limited steps: %.1f%%\n', ...
                tag, SEC_fc_mean, limited_pct);
    else
        fprintf('%s FC-yield: (FC never on)\n', tag);
    end
end


% ---- Fig.4-style polarization (M3) + snapshot ----
i_pol = linspace(0, 3, 301);                 % was 0..2 → now 0..3 A/cm^2
V_pol = arrayfun(@(ii) cell_voltage(ii, pm3), i_pol);
figure('Name','Fig4_Polarization'); 
plot(i_pol, V_pol, 'LineWidth', 2); grid on
xlabel('Current density i [A/cm^2]'); ylabel('Cell voltage V [V]');
title(sprintf('M3 Polarization  @ T=%.0f^\\circC, p_{H2}=p_{O2}=%.2g bar', pm3.el.T-273.15, pm3.el.P_H2));
xlim([0 3]);                                % optional: show full 0..3 range
i_chk = [0, 0.2, 0.5, 1.0, 2.0, 3.0];       % snapshot points extended
V_chk = arrayfun(@(ii) cell_voltage(ii, pm3), i_chk);
fprintf('\nPOLARIZATION SNAPSHOT: V(0)=%.3f | V(0.2)=%.3f | V(0.5)=%.3f | V(1.0)=%.3f | V(2.0)=%.3f | V(3.0)=%.3f\n', V_chk);


% ---- Fig.7 (activation parts) ----
[Vn,Vo,Va_an,Va_ca,~] = cell_voltage_parts(i_pol, pm3);
figure('Name','Fig7_Activation'); 
plot(i_pol,Va_an,'m','LineWidth',2); hold on
plot(i_pol,Va_ca,'b','LineWidth',2); grid on
xlabel('Cell current density [A/cm^2]'); ylabel('Cell voltage [V]');
title('Anode and cathode activation overpotentials (T=80 ^\circC)');
legend('Anode','Cathode','Location','southeast');

% ---- Fig.8 (i0_an sweep) ----
i0_list = [5e-12,1e-11,1e-10,1e-9];
figure('Name','Fig8_i0an'); hold on; grid on
for k=1:numel(i0_list)
    q = pm3; q.el.i0_ref_an = i0_list(k);
    [~,~,~,~, V] = cell_voltage_parts(i_pol, q);
    plot(i_pol, V, 'LineWidth', 2);
end
xlabel('Cell current density [A/cm^2]'); ylabel('Cell voltage [V]');
title('Impact of i0_{an,ref} on polarization (T=80 ^\circC)');
legend(arrayfun(@(v)sprintf('i0_{an,ref}=%.0e A/cm^2',v), i0_list,'uni',0),'Location','southeast');

i_eta = [0, linspace(1e-6, 2.0, 400)];     % include 0 explicitly
eta   = arrayfun(@(ii) faraday_eff(ii, pm3), i_eta);

figure('Name','etaF_vs_i_single');
plot(i_eta, eta, 'LineWidth', 2); grid on
xlabel('Current density i  [A/cm^2]');
ylabel('\eta_F  [-]');
title(sprintf('\\eta_F vs i  @ T=%.0f^\\circC, p_{H2}=%.3g bar, p_{O2}=%.3g bar', ...
      pm3.el.T-273.15, pm3.el.P_H2, pm3.el.P_O2));
ylim([0 1.02]);


% ----------------- helpers -----------------
function [Vnern, Vohm, Vact_an, Vact_ca, Vtot] = cell_voltage_parts(i, p)
  i=i(:)'; R=p.R; F=p.F; T=p.el.T;
  E0 = p.el.DG0_R/(p.z*F);
  Vnern = E0 + (R*T/(p.z*F))*log((p.el.P_H2/1.01325).*sqrt(p.el.P_O2/1.01325)./p.el.a_H2O);
  Vnern = Vnern.*ones(size(i));
  sigma = (F^2*p.el.C_Hp_mol_cm3*p.el.D_Hp_cm2_s)/(R*T);
  Vohm  = (p.el.delta_mem_cm/max(sigma,1e-12)).*i;
  g_an = p.el.phiI_an*p.el.mM_an*(6/(p.el.rhoM_an*p.el.dM_an));
  g_ca = p.el.phiI_ca*p.el.mM_ca*(6/(p.el.rhoM_ca*p.el.dM_ca));
  i0s_an = p.el.i0_ref_an*exp(-p.el.Ea_an/R*(1/T-1/p.el.Tref));
  i0s_ca = p.el.i0_ref_ca*exp(-p.el.Ea_ca/R*(1/T-1/p.el.Tref));
  i0_an = g_an*i0s_an; i0_ca = g_ca*i0s_ca;
  Vact_an = (R*T/(p.el.alpha_an*F))*asinh(i./(2*max(i0_an,1e-20)));
  Vact_ca = (R*T/(p.el.alpha_ca*F))*asinh(i./(2*max(i0_ca,1e-20)));
  Vtot = Vnern + Vohm + Vact_an + Vact_ca;
end
