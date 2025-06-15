
% HESS Model 1 – Simplified Baseline
% Author: Aly Khedr

clear; clc; close all;

%% --- Input Data ---
dt = 1;                      % Time step [hr]
T = 48;                      % Simulation time [hr]
time = 0:dt:T;

fprintf('%5s %10s %10s %15s %15s %15s\n', ...
    'Time', '  P_{el}', 'P_{fc}', ...
    '   m\_dot\_el', '   m\_dot\_fc', '      Tank H_2 [kg]');

%% --- Main Loop ---
for t = 2:length(time)
    % Electrolyzer hydrogen production [kg/h]
    if P_el(t) > 0
        m_dot_el(t) = eta_el * P_el(t) / LHV_H2;
    end

    % Fuel cell hydrogen consumption [kg/h]
    if P_fc(t) > 0
        m_dot_fc(t) = P_fc(t) / (eta_fc * LHV_H2);
    end

    % Update tank level [kg]
    mass_H2(t) = mass_H2(t-1) + (m_dot_el(t) - m_dot_fc(t)) * dt;

    % Enforce tank capacity limits
    mass_H2(t) = min(max(mass_H2(t), 0), max_mass_H2);

    fprintf('%5d %10.2f %10.2f %15.4f %15.4f %15.4f\n', ...
            t-1, P_el(t), P_fc(t), m_dot_el(t), m_dot_fc(t), mass_H2(t));
    
end

%% --- Plot ---
figure;
subplot(2,1,1);
plot(time, P_el, 'b', time, P_fc, 'r');
title('Power Profile'); ylabel('Power [kW]'); legend('P_{el}', 'P_{FC}'); grid on;

subplot(2,1,2);
plot(time, mass_H2, 'k');
ylabel('Tank Level [kg]'); xlabel('Time [hr]');
title('Hydrogen Storage Tank Level'); grid on;
