%% PINN for HESS Level-3 
% Two small PINNs with normalized, unit-consistent physics losses
clear; clc; close all;

%% 1) Load & prepare data
trainTbl = readtable('model3_train_data.csv');
testTbl  = readtable('model3_test_data.csv');
% Derive next-hour tank mass
trainTbl.m_H2_next = [trainTbl.m_H2(2:end); NaN];
testTbl.m_H2_next  = [testTbl.m_H2(2:end); NaN];
% Drop first/last incomplete rows
trainTbl = trainTbl(2:end-1,:);
testTbl  = testTbl(2:end-1,:);

%% 2) Extract inputs and targets, then normalize
Xtr = [trainTbl.P_el, trainTbl.P_fc, trainTbl.m_H2_prev, trainTbl.i_sol, trainTbl.V_cell]';
Xts = [testTbl.P_el,  testTbl.P_fc,  testTbl.m_H2_prev,  testTbl.i_sol,  testTbl.V_cell ]';

Ytr_m = trainTbl.m_H2_next';
Yts_m = testTbl.m_H2_next';

Ytr_f = [ trainTbl.m_dot_el,  trainTbl.m_dot_fc,  trainTbl.eta_el ]';
Yts_f = [  testTbl.m_dot_el,   testTbl.m_dot_fc,   testTbl.eta_el ]';

[xTr,   xPS]  = mapminmax(Xtr);
[yMtr, yMPS]  = mapminmax(Ytr_m);
[yFtr, yFPS]  = mapminmax(Ytr_f);
xTs_n        = mapminmax('apply', Xts, xPS);

%% 3) Define small PINN architectures
% Mass-PINN: 5 -> 8 -> 1
layersM = [
    featureInputLayer(5,'Normalization','none','Name','in')
    fullyConnectedLayer(8,'Name','fc1_m')
    tanhLayer('Name','tanh1_m')
    fullyConnectedLayer(1,'Name','out_m')
];
netM = dlnetwork(layerGraph(layersM));

% Flow-PINN: 5 -> 8 -> 4
layersF = [
    featureInputLayer(5,'Normalization','none','Name','in')
    fullyConnectedLayer(8,'Name','fc1_f')
    tanhLayer('Name','tanh1_f')
    fullyConnectedLayer(3,'Name','out_f')
];
netF = dlnetwork(layerGraph(layersF));

%% 4) Hyperparameters & constants
numEpochs     = 500;
miniBatchSize = min(32, size(xTr,2));
learnRate     = 1e-3;
weightMass    = 1e-3;   % mass-balance loss weight
weightEnergy  = 1e-3;   % energy-balance loss weight

% Convert LHV to kWh/kg for consistent power units
LHV_kWh_per_kg = 120e6 / 3.6e6;    % ≈33.33 kWh/kg

% Reference scales for normalization
P_ref = max(trainTbl.P_el + trainTbl.P_comp);   % kW
M_ref = max(trainTbl.m_H2_next);               % kg

gAvgM = []; gAvgSqM = [];
gAvgF = []; gAvgSqF = [];

%% 5) Custom training loop
numIt = max(floor(size(xTr,2)/miniBatchSize),1);
for epoch = 1:numEpochs
    idx = randperm(size(xTr,2));
    epochLoss = 0;
    for i = 1:numIt
        batchIdx = idx((i-1)*miniBatchSize+1 : min(i*miniBatchSize,end));
        Xb  = dlarray(xTr(:,batchIdx),'CB');
        yMb = dlarray(yMtr(:,batchIdx),'CB');
        yFb = dlarray(yFtr(:,batchIdx),'CB');

        % Compute gradients and loss
        [gM, gF, lossVal] = dlfeval(@computeGradients, netM, netF, ...
            Xb, yMb, yFb, xPS, yMPS, yFPS, LHV_kWh_per_kg, P_ref, M_ref, weightMass, weightEnergy);

        % Update parameters using Adam
        [netM, gAvgM, gAvgSqM] = adamupdate(netM, gM, gAvgM, gAvgSqM, epoch, learnRate);
        [netF, gAvgF, gAvgSqF] = adamupdate(netF, gF, gAvgF, gAvgSqF, epoch, learnRate);

        epochLoss = epochLoss + double(gather(extractdata(lossVal)));
    end
    fprintf('Epoch %3d/%d — avg loss = %.4f\n', epoch, numEpochs, epochLoss/numIt);
end

%% 6) Test set evaluation
Xts_dl  = dlarray(xTs_n,'CB');
yMhat_n = predict(netM, Xts_dl);
yFhat_n = predict(netF, Xts_dl);
Y_mhat  = mapminmax('reverse', extractdata(yMhat_n), yMPS)';
Y_fhat  = mapminmax('reverse', extractdata(yFhat_n), yFPS)';

% True values
Y_mtrue = Yts_m';
Y_ftrue = Yts_f';

% Compute RMSE and R^2
rmse_m = sqrt(mean((Y_mhat - Y_mtrue).^2));
ssRes  = sum((Y_mtrue - Y_mhat).^2);
ssTot  = sum((Y_mtrue - mean(Y_mtrue)).^2);
r2_m   = 1 - ssRes/ssTot;

nFlow = size(Y_ftrue,2);
rmse_f = zeros(nFlow,1);
r2_f   = zeros(nFlow,1);
flowNames = {'m_dot_el','m_dot_fc','eta_el'};
for j = 1:nFlow
    yP = Y_fhat(:,j); yT = Y_ftrue(:,j);
    rmse_f(j) = sqrt(mean((yP - yT).^2));
    r2_f(j)   = 1 - sum((yT - yP).^2)/sum((yT - mean(yT)).^2);
end

% Display results
fprintf('\n--- PINN Performance ---\n');
fprintf('Mass-PINN    : RMSE = %.3f kg, R^2 = %.4f\n', rmse_m, r2_m);
for j = 1:nFlow
    fprintf('Flow-PINN %s: RMSE = %.3f, R^2 = %.4f\n', flowNames{j}, rmse_f(j), r2_f(j));
end

%% 7) Parity plots
allNames = [{'m_H2_next'}, flowNames];
allTrue  = [Y_mtrue, Y_ftrue];
allPred  = [Y_mhat,   Y_fhat];
figure('Position',[100 100 1200 600]);
for k = 1:numel(allNames)
    subplot(2,3,k);
    scatter(allTrue(:,k), allPred(:,k), 20, 'filled'); hold on;
    mn = min([allTrue(:,k); allPred(:,k)]);
    mx = max([allTrue(:,k); allPred(:,k)]);
    plot([mn mx],[mn mx],'k--','LineWidth',1);
    axis([mn mx mn mx]); axis square;
    xlabel('True'); ylabel('Predicted');
    title(allNames{k},'Interpreter','none'); grid on;
end

%% 8) Helper: compute gradients and physics loss
function [gradM, gradF, loss] = computeGradients(...
    netM, netF, X, yM, yF, xPS, yMPS, yFPS, LHV_kWh, P_ref, M_ref, wM, wE)

    % Forward
    pM = forward(netM, X);   % 1×B
    pF = forward(netF, X);   % 4×B

    % Data MSE
    LdM = mse(yM, pM);
    LdF = mse(yF, pF);

    % Un-normalize
    uM = yMPS.gain .* extractdata(pM) + yMPS.xoffset;    % 1×B
    uF = yFPS.gain .* extractdata(pF) + yFPS.xoffset;    % 4×B
    Xu = xPS.gain  .* extractdata(X) + xPS.xoffset;      % 5×B

    % Mass residual (normalized by M_ref)
    prev = Xu(3,:);
    Rm   = (uM - (prev + (uF(1,:)-uF(2,:))*1)) ./ (M_ref + eps);
    Lm   = mean(Rm.^2);

    % Energy residual (all in kW)
    P_el  = Xu(1,:);
    P_chem= uF(1,:) * LHV_kWh;    % kg/h * kWh/kg = kW
    P_fc   = uF(3,:);                      % fuel‑cell power out
    R_E    = (P_el - P_chem + P_fc) ./ (P_ref + eps);
    Le    = mean(R_E.^2);

    % Composite loss
    loss = LdM + LdF + wM*Lm + wE*Le;

    % Gradients
    gradM = dlgradient(loss, netM.Learnables);
    gradF = dlgradient(loss, netF.Learnables);
end
