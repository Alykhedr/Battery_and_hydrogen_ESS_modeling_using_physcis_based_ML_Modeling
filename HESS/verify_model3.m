% verify_model3.m — minimal polarization check (Model-3 only)
close all; clc;

% 1) Build your params with electrochem ON
time = (0:1:1)'; 
p    = params(time, 3);

% 2) Digitized reference points (given)
i_ref   = [0    0.20  0.40  0.60  0.80  1.00]'; 
V25_ref = [1.23 2.30  2.50  2.65  3.85  3.60]'; 
V80_ref = [1.23 1.80  2.05  2.20  2.35  2.80]';

% 3) Two temps using your same params (override only T)
p25 = p; p25.el.T = 298;   % 25°C
p80 = p; p80.el.T = 353;   % 80°C

% 4) Model voltages at the exact ref points
V25_mod = arrayfun(@(x) cell_voltage(x, p25), i_ref);
V80_mod = arrayfun(@(x) cell_voltage(x, p80), i_ref);

% 5) Metrics (RMSE/MAE/Bias) printed to console
e25 = V25_mod - V25_ref;  
e80 = V80_mod - V80_ref;

m25_RMSE = sqrt(mean(e25.^2, 'omitnan'));
m25_MAE  = mean(abs(e25),   'omitnan');
m25_Bias = mean(e25,        'omitnan');

m80_RMSE = sqrt(mean(e80.^2, 'omitnan'));
m80_MAE  = mean(abs(e80),    'omitnan');
m80_Bias = mean(e80,         'omitnan');

fprintf('25°C: RMSE=%.3f V  MAE=%.3f V  Bias=%+.3f V\n', m25_RMSE, m25_MAE, m25_Bias);
fprintf('80°C: RMSE=%.3f V  MAE=%.3f V  Bias=%+.3f V\n', m80_RMSE, m80_MAE, m80_Bias);

% 6) Single plot: model vs reference
figure('Color','w'); hold on; grid on;
plot(i_ref, V25_mod,'-o','LineWidth',1.6,'DisplayName','Model 25°C');
plot(i_ref, V80_mod,'-s','LineWidth',1.6,'DisplayName','Model 80°C');
plot(i_ref, V25_ref,'--^','LineWidth',1.2,'DisplayName','Ref 25°C');
plot(i_ref, V80_ref,'--v','LineWidth',1.2,'DisplayName','Ref 80°C');
xlabel('Current density i (A/cm^2)'); ylabel('Cell voltage V (V)');
title('PEM Electrolyzer Polarization: Model vs Reference');
legend('Location','northwest');

