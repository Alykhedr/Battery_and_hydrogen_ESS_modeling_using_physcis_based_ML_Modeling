%% One-Net FFNN: predict [m_dot_el(t), m_H2_next] from [P_el, P_fc, m_H2, m_H2_prev]
clear; clc; close all;

%% 0) RNG
seed = 42;
rng(seed,'twister');

%% 1) Load
trainTbl = readtable('model3_train_data.csv');
testTbl  = readtable('model3_test_data.csv');

%% 2) Build next-step label & cleanup
trainTbl.m_H2_next = [trainTbl.m_H2(2:end); NaN];
testTbl.m_H2_next  = [testTbl.m_H2(2:end);  NaN];

% Drop the last row (NaN next) and any accidental NaNs in inputs
keepVars = {'Hour','P_el','P_fc','m_H2','m_H2_prev','m_dot_el','m_H2_next'};
trainTbl = rmmissing(trainTbl(:, keepVars));
testTbl  = rmmissing(testTbl(:,  keepVars));

%% 3) Define inputs/targets
X_cols = {'P_el','P_fc','m_H2','m_H2_prev'};
Y_cols = {'m_dot_el','m_H2_next'};

X_tr = trainTbl{:, X_cols}';   % [nFeat x nTrain]
Y_tr = trainTbl{:, Y_cols}';   % [2 x nTrain]
X_ts = testTbl{:,  X_cols}';
Y_ts = testTbl{:,  Y_cols}';

%% 4) Normalise (fit on train only)
[Xn_tr, psX] = mapminmax(X_tr);           % feature-wise
[Yn_tr, psY] = mapminmax(Y_tr);           % target-wise
Xn_ts = mapminmax('apply', X_ts, psX);    % apply train scaling

%% 5) Train/Val split (90/10)
% Sort by time if you have Hour or a timestamp; else use row order
n = size(Xn_tr,2);
cut = round(0.9*n);
trainInd = 1:cut;
valInd   = cut+1:n;


%% 6) Network
% A small, low-complexity net; adjust widths if needed
net = feedforwardnet([20,10], 'trainscg');
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.valInd   = valInd;
net.divideParam.testInd  = [];

% Regularization & training display
net.performParam.regularization = 0.1;
net.trainParam.showWindow = false;
net.trainParam.max_fail   = 12;   % early stopping patience

% Element-wise loss weights (2 x nTrain) to emphasise targets
% Give higher weight to m_H2_next; still weight m_dot_el > 1
w_mdot   = 1.7;
w_mH2nxt = 2.5;
EW = ones(2, n);                 % default 1
EW(1, :) = w_mdot;               % row 1 -> m_dot_el
EW(2, :) = w_mH2nxt;             % row 2 -> m_H2_next

%% 7) Train
[net, tr] = train(net, Xn_tr, Yn_tr, [], [], EW);

% Plot convergence
figure('Name','One-Net Convergence');
plot(tr.perf,'b-'); hold on; plot(tr.vperf,'r--'); grid on;
xlabel('Epoch'); ylabel('MSE (weighted)'); legend('Train','Val');
title('One-Net: Training vs Validation');



%% 8) Predict (test) & invert scaling
Yn_pred_ts = net(Xn_ts);
Y_pred_ts  = mapminmax('reverse', Yn_pred_ts, psY)';   % [nTest x 2]
Y_true_ts  = Y_ts';

% Split back
y_mdot_true = Y_true_ts(:,1);  y_mdot_pred = Y_pred_ts(:,1);
y_mH2n_true = Y_true_ts(:,2);  y_mH2n_pred = Y_pred_ts(:,2);

% Split back into train/val using the same indices you used to train
Xn_tr_all = Xn_tr; Yn_tr_all = Yn_tr;
Xn_tr_sub = Xn_tr_all(:, trainInd);
Xn_va_sub = Xn_tr_all(:, valInd);

% Predictions
Yp_tr  = mapminmax('reverse', net(Xn_tr_sub), psY)';
Yp_va  = mapminmax('reverse', net(Xn_va_sub), psY)';
Yp_ts  = mapminmax('reverse', net(Xn_ts),    psY)';

Yt_tr  = mapminmax('reverse', Yn_tr_all(:,trainInd), psY)';  % de-normalize truths
Yt_va  = mapminmax('reverse', Yn_tr_all(:,valInd),   psY)';

rmse  = @(a,b) sqrt(mean((a-b).^2));
r2    = @(a,b) 1 - sum((a-b).^2)/sum((a-mean(a)).^2);

names = {'m\_dot\_el','m\_H2\_next'};

for j = 1:2
    fprintf('\nTarget: %s\n', names{j});
    fprintf('  TRAIN: RMSE=%.4g  R2=%.6f\n', rmse(Yt_tr(:,j),Yp_tr(:,j)), r2(Yt_tr(:,j),Yp_tr(:,j)));
    fprintf('  VAL  : RMSE=%.4g  R2=%.6f\n', rmse(Yt_va(:,j),Yp_va(:,j)), r2(Yt_va(:,j),Yp_va(:,j)));
    fprintf('  TEST : RMSE=%.4g  R2=%.6f\n', rmse(Y_true_ts(:,j),Y_pred_ts(:,j)), r2(Y_true_ts(:,j),Y_pred_ts(:,j)));
end


%% 9) Metrics
rmse  = @(a,b) sqrt(mean((a-b).^2));
r2    = @(a,b) 1 - sum((a-b).^2)/sum((a-mean(a)).^2);
mae   = @(a,b) mean(abs(a-b));

stats = table( ...
    ["m_dot_el"; "m_H2_next"], ...
    [rmse(y_mdot_true,y_mdot_pred);  rmse(y_mH2n_true,y_mH2n_pred)], ...
    [mae(y_mdot_true,y_mdot_pred);   mae(y_mH2n_true,y_mH2n_pred)], ...
    [r2(y_mdot_true,y_mdot_pred);    r2(y_mH2n_true,y_mH2n_pred)], ...
    'VariableNames', {'Target','RMSE','MAE','R2'});

fprintf('\n--- One-Net Test Performance ---\n'); disp(stats);

%% 10) Parity & time series
figure('Name','Parity','Position',[100 100 900 400]);
subplot(1,2,1);
scatter(y_mdot_true, y_mdot_pred, 18, 'filled'); hold on;
lm = [min([y_mdot_true; y_mdot_pred]) max([y_mdot_true; y_mdot_pred])];
plot(lm,lm,'k--'); axis square; grid on;
xlabel('True m\_dot\_el'); ylabel('Pred'); title('m\_dot\_el');

subplot(1,2,2);
scatter(y_mH2n_true, y_mH2n_pred, 18, 'filled'); hold on;
lm = [min([y_mH2n_true; y_mH2n_pred]) max([y_mH2n_true; y_mH2n_pred])];
plot(lm,lm,'k--'); axis square; grid on;
xlabel('True m\_H2\_next'); ylabel('Pred'); title('m\_H2\_next');

figure('Name','Time Series','Position',[100 100 1000 500]);
t = testTbl.Hour;
subplot(2,1,1);
plot(t, y_mdot_true,'b-', t, y_mdot_pred,'r--','LineWidth',1.3); grid on;
ylabel('m\_dot\_el'); legend('True','Pred','Location','best');
title('Electrolyser H2 Production');

subplot(2,1,2);
plot(t, y_mH2n_true,'b-', t, y_mH2n_pred,'r--','LineWidth',1.3); grid on;
xlabel('Hour'); ylabel('m\_H2\_next'); legend('True','Pred','Location','best');
title('Tank H2 Mass (next step)');
