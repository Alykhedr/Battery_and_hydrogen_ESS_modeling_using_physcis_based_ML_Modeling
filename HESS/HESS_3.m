%  3rd level HESS Model 3
% Author: Aly Khedr 
% Description: 

clear; clc; close all;

%% ------------------ Input Parameters ------------------
dt = 1;                      % Time step [hr]
T = 48;                     % Total simulation time [hr]
time = 0:dt:T;

% Universal Constants
F = 96485;                  % Faraday constant [C/mol]
R = 8.314;                  % Gas constant [J/mol/K]
T_el = 350;                 % Electrolyzer temperature [K]
T_tank = 300;               % Tank temperature [K]
LHV_H2 = 33.33;             % Lower heating value [kWh/kg]
M_H2 = 0.002016;            % Hydrogen molar mass [kg/mol]

% Electrolyzer Parameters
V_el_r = 1.23;              % Reversible voltage [V]
Omega_el = 0.2;             % Ohmic resistance [Ohm·cm^2]
sigma_p = 0.5; sigma_n = 0.5;
ip = 1e-4; in_ = 1e-3;       % Exchange current densities [A/cm^2]
A_el = 1000;                % Active area [cm^2]
N_el = 400;                 % Number of cells
L = 0.0178;                 % Membrane thickness [cm]

% Permeability
epsilon_H2f = 2e-10;
epsilon_O2f = 1e-10;
epsilon_H2p = 1e-10;
xi_H2 = 0.01; xi_O2 = 0.01;
P_H2 = 5; P_O2 = 5;          % [bar]

% Fuel Cell Parameters
eta_FC = 0.60;

% Compressor Parameters
Cp_H2 = 14.3; eta_C = 0.75; T_C = 300; gamma = 1.4;
P_inlet = 10e5; P_outlet = 70e5;

% Tank Parameters
V_tank = 12;                % [m^3]

%% ------------------ Initialization ------------------
n_H2 = zeros(1, length(time));
P_el = zeros(1, length(time));
P_fc = zeros(1, length(time));
m_dot_el = zeros(1, length(time));
m_dot_fc = zeros(1, length(time));
P_tank = zeros(1, length(time));
P_C = zeros(1, length(time));
i_sol = zeros(1, length(time));
V_cell = zeros(1, length(time));
eta_tf = zeros(1, length(time));

n_H2(1) = 25 / M_H2;  % initial H2 [mol]

%% ------------------ Power Profile ------------------
max_power = 400;
rng(1);
for t = 2:length(time)
    if rand < 0.5
        P_el(t) = rand * max_power;
    else
        P_fc(t) = rand * max_power;
    end
end

%% ------------------ Simulation Loop ------------------
fprintf('%5s %10s %10s %15s %15s %15s\n', 'Time', 'P_el', 'P_fc', 'H2 Mass [kg]', 'P_C [kW]', 'P_tank [bar]');

for t = 2:length(time)
    Pel_t = P_el(t);
    Pfc_t = P_fc(t);

    %% Electrolyzer
    if Pel_t > 0
        % Solve for current density using fzero
        fun = @(i) N_el * i * A_el / 1000 * ...
            (V_el_r + (R*T_el/(2*F))*log(P_H2*sqrt(P_O2)) + Omega_el*i + ...
            (R*T_el)/(sigma_p*F)*asinh(i/(2*ip)) + (R*T_el)/(sigma_n*F)*asinh(i/(2*in_))) - Pel_t;
        i_guess = fzero(fun, [0.01 2]);
        i_sol(t) = i_guess;

        % Voltage
        V_cell(t) = V_el_r + (R*T_el/(2*F))*log(P_H2*sqrt(P_O2)) + Omega_el*i_guess + ...
                    (R*T_el)/(sigma_p*F)*asinh(i_guess/(2*ip)) + ...
                    (R*T_el)/(sigma_n*F)*asinh(i_guess/(2*in_));

        % Pressures
        pi_H2n = P_H2 + xi_H2 * (1 + epsilon_H2p / epsilon_H2f) * i_guess;
        pi_O2p = P_O2 + xi_O2 * i_guess;

        % Faraday Efficiency
        eta = 1 - (2 * F * epsilon_H2f * pi_H2n * L + ...
                  4 * F * epsilon_O2f * pi_O2p * L + ...
                  2 * F * epsilon_H2p * pi_H2n * pi_O2p * L);
        eta_tf(t) = max(min(eta, 1), 0);

        % Hydrogen Production
        m_dot_el(t) = N_el * eta_tf(t) * i_guess * A_el / (2*F) * M_H2 * 3600;
    end

    %% Fuel Cell
    if Pfc_t > 0
        required_mass = Pfc_t / (eta_FC * LHV_H2);
        available_mass = n_H2(t-1) * M_H2;
        m_dot_fc(t) = min(required_mass, available_mass / dt);
    end

    %% Tank State
    n_dot_el = m_dot_el(t) / M_H2;
    n_dot_fc = m_dot_fc(t) / M_H2;
    n_H2(t) = n_H2(t-1) + (n_dot_el - n_dot_fc) * dt;
    n_H2(t) = max(min(n_H2(t), 50 / M_H2), 0);
    P_tank(t) = (n_H2(t) * R * T_tank) / V_tank;

    %---- 4) Compressor Work ----%
   n_dot_el = m_dot_el(t) / M_H2 / 3600;  

    if n_dot_el > 0 && P_tank(t) < P_outlet
        PR = P_outlet / max(P_tank(t), P_inlet);
        P_C(t) = n_dot_el * Cp_H2 * T_C / eta_C * (PR^((gamma-1)/gamma) - 1);
    else
       P_C(t) = 0;
    end

    fprintf('%5d %10.2f %10.2f %15.4f %15.4f %15.2f\n', ...
        time(t), Pel_t, Pfc_t, n_H2(t)*M_H2, P_C(t)/1000, P_tank(t)/1e5);
end

%% ------------------ Post-Processing ------------------
mass_H2 = n_H2 * M_H2;
P_C_kW = P_C / 1000;
P_tank_bar = P_tank / 1e5;

%% ------------------ Visualization (Model 3, trimmed) ------------------
figure('Position',[100 100 800 600]);

% 1) Power Profile
subplot(2,2,1);
plot(time, P_el, 'b', time, P_fc, 'r', 'LineWidth',1.2);
title('Power Profile'); ylabel('kW'); xlabel('Time [hr]');
legend('P_{el}','P_{fc}'); grid on;

% 2) H2 Mass in Tank
subplot(2,2,2);
plot(time, mass_H2, 'k', 'LineWidth',1.5);
title('H₂ Mass in Tank'); ylabel('kg'); xlabel('Time [hr]');
grid on;

% 3) Compressor Power & Tank Pressure
subplot(2,2,3);
yyaxis left
plot(time, P_C_kW, 'm', 'LineWidth',1.5);
ylabel('Compressor Power [kW]');
yyaxis right
plot(time, P_tank_bar, 'g--', 'LineWidth',1.5);
ylabel('Tank Pressure [bar]');
xlabel('Time [hr]');
title('Compressor & Pressure');
legend('P_C','P_{tank}','Location','best');
grid on;

% 4) Current Density
subplot(2,2,4);
plot(time, i_sol, 'b', 'LineWidth',1.2);
title('Electrolyzer Current Density'); 
ylabel('A/cm^2'); xlabel('Time [hr]');
grid on;
