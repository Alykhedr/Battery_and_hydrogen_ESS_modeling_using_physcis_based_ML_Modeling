% HESS Model 1
% Author: Aly Khedr

clear; clc; close all;

% 1) load profile
[time, P_el, P_fc] = load_profiles();

% 2) grab constants + initial state
p        = params(time,1);
p.enable.compressor = false;
m_dot_el = p.init.m_dot_el;
m_dot_fc = p.init.m_dot_fc;
mass_H2  = p.init.mass_H2;

% 3) compute dt
dt = time(2) - time(1);

% 4) run the HESS loop
for k = 2:numel(time)
  if P_el(k) > 0
    m_dot_el(k) = p.eta_el * P_el(k) / p.LHV_H2;       % kg/h
  end
  if P_fc(k) > 0
    m_dot_fc(k) = P_fc(k) / (p.eta_fc * p.LHV_H2);     % kg/h
  end

  mass_H2(k) = mass_H2(k-1) + (m_dot_el(k-1) - m_dot_fc(k-1))*dt;
  mass_H2(k) = min(max(mass_H2(k),0), p.max_mass_H2);
end

% 5) plot
plot_hess(time, P_el, P_fc, mass_H2, m_dot_el, m_dot_fc);

%% ============================================================
%  Verification vs vendor spec (SYSTEM basis, LHV)
%  - Model-1 has no compressor, so Pin_system = P_el (kW)
%  - SEC_model = Pin_system / m_dot_el  [kWh/kg]
%  - Reference: SEC_ref = 53 kWh/kg (adjust if you prefer 50)
% =============================================================
SEC_ref = 55.54;                           % <-- adjust here if needed

% Make sure model outputs are column vectors, matching time length
P_system_kW = max(0, P_el(:));            % system input power (kW)
m_dot_el    = m_dot_el(:);                % kg/h (from loop)
SEC_mod     = P_system_kW ./ max(m_dot_el, 1e-12);  % kWh/kg

% Only compare above a reasonable load (e.g., >=30% of nominal)
if isfield(p,'P_el_rated') && ~isempty(p.P_el_rated)
    P_nom_kW = p.P_el_rated;
else
    P_nom_kW = max(P_system_kW);
end
min_load_frac = 0.30;
valid = (m_dot_el > 1e-6) & (P_system_kW >= min_load_frac*P_nom_kW);

% Metrics
err  = SEC_mod(valid) - SEC_ref;
MAE  = mean(abs(err));
RMSE = sqrt(mean(err.^2));
BIAS = mean(err);
MAPE = mean(abs(err)/SEC_ref)*100;

% Implied efficiencies (LHV basis)
eta_ref   = p.LHV_H2 / SEC_ref;           % fraction
eta_mod   = p.LHV_H2 ./ SEC_mod(valid);
eta_mean  = mean(eta_mod);
eta_med   = median(eta_mod);

fprintf('\nSEC Verification (Model-1, SYSTEM/LHV)\n');
fprintf('  Points used: %d / %d  (>= %.0f%% load)\n', nnz(valid), numel(valid), 100*min_load_frac);
fprintf('  SEC_ref = %.2f kWh/kg | LHV = %.2f kWh/kg -> eta_ref = %.1f%%\n', SEC_ref, p.LHV_H2, 100*eta_ref);
fprintf('  MAE = %.2f kWh/kg     Bias = %+0.2f\n', MAE, BIAS);


% Safe sample print (use time if aligned, otherwise use index)
t_disp = (1:numel(P_system_kW))';
if numel(time) == numel(P_system_kW)
    t_disp = time;
end

idx = find(valid);
idx = idx(1:min(6, numel(idx)));  % first few valid samples
if ~isempty(idx)
    fprintf('\n   t | Pin_sys | mdot | SEC_mod |  ΔSEC\n');
    fprintf('-----+---------+------+---------+------\n');
    for j = idx(:)'
        fprintf('%4.0f | %7.1f | %4.2f | %7.2f | %+5.2f\n', ...
            t_disp(j), P_system_kW(j), m_dot_el(j), SEC_mod(j), SEC_mod(j)-SEC_ref);
    end
end

