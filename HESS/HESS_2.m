
% 2nd Level HESS Model
% Author: Aly Khedr
% Description: HESS model with electrolyzer, fuel cell, compressor, and ideal gas law-based tank

clear; clc; close all;

%% ------------------ Input Parameters ------------------
% Time settings
dt = 1;                      % Time step [hr]
T = 48;                      % Total time [hr]
time = 0:dt:T;

% Electrolyzer and Fuel Cell Parameters
LHV_H2   = 33.33;            % Lower heating value of H2 [kWh/kg]
eta_el   = 0.75;             % Electrolyzer efficiency
eta_FC   = 0.60;             % Fuel cell efficiency

% Compressor parameters
Cp_H2    = 14.3;             % Specific heat [J/mol·K]
eta_C    = 0.75;             % Compressor efficiency
T_C      = 300;              % Compressor temperature [K]
gamma    = 1.4;              % Adiabatic index

% Ideal gas law constants
R        = 8.314;            % Universal gas constant [J/mol·K]
T_tank   = 300;              % Tank temperature [K]
V_tank   = 12;              % Tank volume [m³]
M_H2     = 0.002016;         % H2 molar mass [kg/mol]

% Compressor pressure range
P_inlet  = 10e5;              % Inlet pressure [Pa]
P_outlet = 70e5;             % Outlet pressure [Pa]

%% ------------------ Initialization ------------------
n_H2      = zeros(1, length(time));     % H2 content [mol]
P_tank    = zeros(1, length(time));     % Tank pressure [Pa]
P_C       = zeros(1, length(time));     % Compressor power [W]
P_el      = zeros(1, length(time));     % Electrolyzer power [kW]
P_fc      = zeros(1, length(time));     % Fuel cell power [kW]
m_dot_el  = zeros(1, length(time));     % H2 production [kg/h]
m_dot_fc  = zeros(1, length(time));     % H2 consumption [kg/h]

% Initial hydrogen state
initial_mass = 0.5 * 50;               % 60% of 200 kg
n_H2(1) = initial_mass / M_H2;

%% ------------------ Power Profile Generation ------------------
max_power = 400;  % kW
rng(1);           % Seed for reproducibility

for t = 2:length(time)
    if rand < 0.5
        P_el(t) = rand * max_power;
        P_fc(t) = 0;
    else
        P_el(t) = 0;
        P_fc(t) = rand * max_power;
    end
end

%% ------------------ Simulation Loop ------------------
fprintf('%5s %10s %10s %15s %15s %15s\n', 'Time', 'P_el [kW]', 'P_fc [kW]', 'Tank H2 [kg]', 'P_C [kW]',' Tank pressure[bar]');

for t = 2:length(time)
    Pel_t = P_el(t);
    Pfc_t = P_fc(t);

    % H2 mass flow rates [kg/h]
    if Pel_t > 0
        m_dot_el(t) = (Pel_t * eta_el) / LHV_H2;
    end
    if Pfc_t > 0
        m_dot_fc(t) = Pfc_t / (eta_FC * LHV_H2);
    end

    % Convert to mol/h
    n_dot_el = m_dot_el(t) / M_H2;
    n_dot_fc = m_dot_fc(t) / M_H2;

    % Update hydrogen content
    n_H2(t) = n_H2(t-1) + (n_dot_el - n_dot_fc) * dt;
    n_H2(t) = min(max(n_H2(t), 0), 50 / M_H2);  % Limit to 0–50 kg range

    % Tank pressure [Pa]
    P_tank(t) = (n_H2(t) * R * T_tank) / V_tank;

    %---- 4) Compressor Work ----%
   n_dot_el = m_dot_el(t) / M_H2 / 3600;  

    if n_dot_el > 0 && P_tank(t) < P_outlet
        PR = P_outlet / max(P_tank(t), P_inlet);
        P_C(t) = n_dot_el * Cp_H2 * T_C / eta_C * (PR^((gamma-1)/gamma) - 1);
    else
       P_C(t) = 0;
    end

    % Step-by-step output
    fprintf('%5d %10.2f %10.2f %15.4f %15.4f %15.2f\n', ...
    time(t), Pel_t, Pfc_t, n_H2(t)*M_H2, P_C(t)/1000, P_tank(t)/1e5);

end

%% ------------------ Post-Processing ------------------
mass_H2     = n_H2 * M_H2;         % H2 mass in kg
P_C_kW      = P_C / 1000;          % Compressor power in kW
P_tank_bar  = P_tank / 1e5;        % Pressure in bar

%% ------------------ Visualization ------------------

figure;

% Plot 1: Power Profiles
subplot(3,1,1);
plot(time, P_el, 'b', 'LineWidth', 1.2); hold on;
plot(time, P_fc, 'r', 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Power [kW]');
title('Electrolyzer and Fuel Cell Power'); grid on;
legend('P_{el}', 'P_{FC}');

% Plot 2: Hydrogen Mass in Tank
subplot(3,1,2);
plot(time, mass_H2, 'k', 'LineWidth', 1.5);
xlabel('Time [hr]'); ylabel('Hydrogen Mass [kg]');
title('Hydrogen Tank Mass'); grid on;

% Plot 3: Compressor Power and Tank Pressure
subplot(3,1,3);
yyaxis left
plot(time, P_C_kW, 'm', 'LineWidth', 1.5);
ylabel('Compressor Power [kW]');
yyaxis right
plot(time, P_tank_bar, 'g--', 'LineWidth', 1.5);
ylabel('Tank Pressure [bar]');
xlabel('Time [hr]');
title('Compressor Power and Tank Pressure');
legend('P_C', 'P_{tank}'); grid on;
