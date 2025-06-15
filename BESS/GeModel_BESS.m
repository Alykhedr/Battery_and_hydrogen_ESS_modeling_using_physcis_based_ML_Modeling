%% soc_static_trap_profile_with_fade.m
% GE‐style SOC evolution for fixed γ values, using a trapezoidal C‐rate profile
% Enhanced: plots both normalized SoC and absolute available capacity to
% visualize capacity fade.
clear; clc; close all;

%% 1 ── Cell & model parameters
Q_N        = 1.1;        % Nominal capacity [Ah]
I_N        = 1.0;        % Rated current  [A]
T_nom      = 25;         % Nominal temperature [°C]
T_amb      = 22;         % Ambient temperature [°C]
n          = -0.03;      % Efficiency exponent
SOC0       = 1.0;        % Initial SoC (100% of faded capacity)

%% 2 ── Time grid
dt_hr      = 1;          % Time step [h]
t_end      = 48;         % Total simulation time [h]
T_vec      = 0:dt_hr:t_end;
N          = numel(T_vec);

%% 3 ── Trapezoidal current‐rate profile (in C‐rates)
pattern_C  = [ 1, 0, -0.5, -0.5, 0 ];  % Discharge 1C, rest, charge 0.5C×2h, rest
period_h   = numel(pattern_C);
I_C        = repmat(pattern_C,1,ceil(N/period_h));
I_C        = I_C(1:N);
I          = I_C * I_N;  % Convert C-rate to current [A]

%% 4 ── Fixed γ states from Preger/Ge data
cycle_count = [   0, 1000, 2000, 3000];
gamma_vals  = [1.00, 0.96, 0.94, 0.92];
colors      = lines(numel(cycle_count));

%% 5 ── Simulation & data storage
SOC_matrix       = zeros(numel(gamma_vals), N);
Q_available_mat  = zeros(numel(gamma_vals), N);

for j = 1:numel(cycle_count)
    gamma = gamma_vals(j);
    SOC   = zeros(1,N);
    SOC(1)= SOC0;
    
    for k = 2:N
        % instantaneous temperature factor\        
        K_T = 1 + 0.008 * (T_nom - T_amb);
        I_k = I(k);
        % coulombic efficiency
        if I_k ~= 0
            eta = (abs(I_k)/I_N)^n;
            eta = min(max(eta,0),1);
        else
            eta = 1;
        end
        % SoC update
        dSOC = (K_T * eta * I_k * dt_hr) / (gamma * Q_N);
        SOC(k) = SOC(k-1) - dSOC;
        SOC(k) = min(max(SOC(k),0),1);
    end
    % store normalized SoC
    SOC_matrix(j,:) = SOC;
    % compute absolute available capacity [Ah]
    Q_available_mat(j,:) = gamma * Q_N .* SOC;
end

%% 6 ── Visualization
figure('Position',[100 100 800 600]);

% Subplot 1: normalized SoC
subplot(2,1,1); hold on; grid on;
for j = 1:numel(cycle_count)
    plot(T_vec, SOC_matrix(j,:)*100, 'LineWidth',1.5, 'Color',colors(j,:), ...
         'DisplayName', sprintf('%d cycles (γ=%.2f)', cycle_count(j), gamma_vals(j)));
end
xlabel('Time [h]');
ylabel('SoC [% of faded capacity]');
title('Normalized SoC Evolution under Capacity Fade');
legend('Location','southwest');
set(gca,'FontSize',12);

% Subplot 2: absolute available capacity
subplot(2,1,2); hold on; grid on;
for j = 1:numel(cycle_count)
    plot(T_vec, Q_available_mat(j,:), 'LineWidth',1.5, 'Color',colors(j,:), ...
         'DisplayName', sprintf('γ=%.2f → %.2f Ah', gamma_vals(j), gamma_vals(j)*Q_N));
end
xlabel('Time [h]');
ylabel('Available Capacity [Ah]');
title('Absolute Available Capacity under Capacity Fade');
legend('Location','southwest');
set(gca,'FontSize',12);
