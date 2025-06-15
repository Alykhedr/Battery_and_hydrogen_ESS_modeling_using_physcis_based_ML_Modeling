 % 2nd level HESS Model
 % Author: Aly Khedr 
 % Description: HESS model wihth detailed electrolyzer electrochmeical
 % model
clear; clc; close all;

%% --- Input Data ---
dt = 1;                      % Time step [hr]
T = 48;                      % Simulation time [hr]
time = 0:dt:T;

% Constants
F = 96485;                   % Faraday constant [C/mol]
R = 8.314;                   % Universal gas constant [J/mol/K]
T_el = 350;                  % Electrolyzer temperature [K]
LHV_H2 = 33.33;              % [kWh/kg]
H_H2_h = 39.4e6;             % High heating value [J/kg]
M_H2 = 0.002016;             % H2 molar mass [kg/mol]
Z = 1.1;                     %compressibility factor

% Electrolyzer parameters
V_el_r = 1.23;               % Reversible voltage [V]
beta_H2O = 1;                % Water activity
Omega_el = 0.2;              % Internal resistance [Ohm·cm^2]
sigma_p = 0.5; sigma_n = 0.5;% Transfer coefficients
ip = 1e-4; in_ = 1e-3;       % Exchange current densities [A/cm^2]
A_el = 1000;                 % Active area [cm^2]
N_el = 400;                  % Number of cells
k_el = 1;                    % Unit scaling

% Gas permeabilities and membrane
epsilon_H2 = 2e-10; epsilon_O2 = 1e-10; epsilon_H2p = 1e-10; % [mol/cm/s/bar]
L = 0.0178;                  % Membrane thickness [cm]
xi_H2 = 0.01; xi_O2 = 0.01;  % Fitting params [bar·cm^2/A]
P_H2 = 5; P_O2 = 5;          % Electrolyzer output pressure [bar]

% Compressor constants
 eta_C = 0.75;gamma = 1.4;
P_inlet = 5e5; P_outlet = 35e5;

% Tank constants
T_tank = 300; V_tank = 7.5; max_mass_H2 = 200; % kg

%% === Initialize Variables ===
n_H2 = zeros(1, length(time));
P_tank = zeros(1, length(time));
P_C = zeros(1, length(time));
P_el = zeros(1, length(time));
P_fc = zeros(1, length(time));
m_dot_el = zeros(1, length(time));
m_dot_fc = zeros(1, length(time));
eta_el = zeros(1,length(time));
eta_F = zeros(1,length(time));

% Initial mass and pressure
n_H2(1) = 0.6 * max_mass_H2 / M_H2;
P_tank(1) = 20; 

%% === Generate Random Power Profile ===
max_power = 800; rng(1);
for t = 2:length(time)
    if rand < 0.5
        P_el(t) = rand * max_power;
        P_fc(t) = 0;
    else
        P_el(t) = 0;
        P_fc(t) = rand * max_power;
    end
end

%% === Main Loop ===
for t = 2:length(time)
    Pel = P_el(t); Pfc = P_fc(t);

    % Fuel Cell
    if Pfc > 0
        m_dot_fc(t) = Pfc / (0.60 * LHV_H2);
    end

    % Electrolyzer model if active
    if Pel > 0
        % Assume constant current density i
        i = 1.0;  % A/cm^2 

        % Pressures
        PH2_n = P_H2 + xi_H2 * (1 + epsilon_H2p / epsilon_H2) * i;
        PO2_p = P_O2 + xi_O2 * i;

        % Voltage calculation
        V_log = R * T_el / (2 * F) * log((P_H2 * sqrt(P_O2)) / beta_H2O);
        V_ohmic = Omega_el * i;
        V_anode = (R * T_el / (sigma_p * F)) * asinh(i / (2 * ip));
        V_cathode = (R * T_el / (sigma_n * F)) * asinh(i / (2 * in_));
        V_el = V_el_r + V_log + V_ohmic + V_anode + V_cathode;

        % Power
        P_el(t) = N_el * V_el * i * A_el / 1000; % [W] to [kW]

        % Faraday efficiency
        eta_F = 1 - (2 * F * epsilon_H2 * PH2_n + 4 * F * epsilon_O2 * PO2_p + ...
            2 * F * epsilon_H2p * (PH2_n - PO2_p)) / L;

        % Hydrogen production
        m_dot_el(t) = ((k_el * N_el * eta_F * i * A_el) / (2 * F))* M_H2 * 3600;  % kg/h
      
        % Electrolyzer efficiency
        eta_el(t) = (m_dot_el(t)/3600 * H_H2_h) / (P_el(t)*1000);  
    end

    % Convert to mol/h
    n_dot_el = m_dot_el(t) / M_H2;
    n_dot_fc = m_dot_fc(t) / M_H2;

    % Update tank content
    n_H2(t) = n_H2(t-1) + (n_dot_el - n_dot_fc);
    n_H2(t) = min(max(n_H2(t), 0), max_mass_H2 / M_H2);

    % Tank pressure
     P_tank(t) = P_tank(t-1) + ...
    (Z * R * T_tank / max_mass_H2) * (m_dot_el(t) - m_dot_fc(t)) * dt;

    % Apply safety limits (20–50 bar)
    P_tank(t) = min(max(P_tank(t), 20), 50);
    % Compressor
     if m_dot_el(t) > 0
        P_C(t) = (m_dot_el(t) * R * T_el) / (eta_C * (gamma - 1)) * ...
                 ((P_outlet / P_inlet)^((gamma - 1)/gamma) - 1) / 1000;  % [kW]
    else
        P_C(t) = 0;
    end


end

%% === Plot ===
P_C_kW = P_C ;
mass_H2 = n_H2 * M_H2;

figure;
subplot(3,1,1);
plot(time, P_el, 'b', time, P_fc, 'r');
title('Power Profile'); ylabel('Power [kW]'); legend('P_{el}', 'P_{FC}'); grid on;

subplot(3,1,2);
plot(time, mass_H2, 'k');
title('Hydrogen Mass in Tank'); ylabel('Mass [kg]'); grid on;

subplot(3,1,3);
yyaxis left; plot(time, P_C_kW, 'm'); ylabel('Compressor Power [kW]');
yyaxis right; plot(time, P_tank, 'g--'); ylabel('Tank Pressure [bar]');
title('Compressor & Tank'); xlabel('Time [hr]'); legend('P_C', 'P_{tank}'); grid on;







