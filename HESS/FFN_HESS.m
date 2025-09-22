%% HESS FFNN surrogate (multi-output, no lags; predict next-step mass + flows)
clear; clc;

% ---------------- IO ----------------
csvPath   = 'model3_train_data.csv';  
modelsDir = 'models_hess';
plotsDir  = 'plots_hess';
if ~exist(modelsDir,'dir'), mkdir(modelsDir); end
if ~exist(plotsDir,'dir'), mkdir(plotsDir); end
folderPath = plotsDir;  % export paths 

% ---------------- Load ----------------
T = readtable(csvPath);
% Expect at least these columns in CSV:
% {'Hour','P_el','P_fc','m_dot_el','m_dot_fc','eta_el','i_model','V','eta_F','P_comp','m_H2','P_tank_bar'}

% Build next-step mass target and trim last row 
T.m_H2_next = [T.m_H2(2:end); T.m_H2(end)];  % temporary
T(end,:) = [];                                % drop last row so next is valid

% ---------------- Features / Targets ----------------
featNames  = {'Hour','P_el','P_fc','m_H2','V','i_model'};   % 6 inputs
targetCols = {'m_H2_next','m_dot_el','m_dot_fc','P_comp','P_tank_bar'}; % 5 outputs

X_tbl = T(:, featNames);
Y_tbl = T(:, targetCols);

% ensure numeric arrays
X_mat = table2array(X_tbl);
Y_mat = table2array(Y_tbl);

% handle NaNs just in case
X_mat(isnan(X_mat)) = 0;
Y_mat(isnan(Y_mat)) = 0;

% ---------------- Metrics helpers ----------------
r2col = @(yPred,yTrue) 1 - sum((yPred - yTrue).^2) ./ max(sum((yTrue - mean(yTrue)).^2), eps);
% aggregate metrics across 5 outputs (mean over columns)
aggRMSE = @(E) sqrt(mean(mean(E.^2,1),2));
aggMAE  = @(E) mean(mean(abs(E),1),2);
aggR2   = @(Yp,Yt) mean(r2col(Yp, Yt), 2);

% ---------------- Training setup ----------------
numRuns = 2;

rmseTrainAll = zeros(numRuns,1);
maeTrainAll  = zeros(numRuns,1);
r2TrainAll   = zeros(numRuns,1);

rmseValAll   = zeros(numRuns,1);
maeValAll    = zeros(numRuns,1);
r2ValAll     = zeros(numRuns,1);

rmseTestAll  = zeros(numRuns,1);
maeTestAll   = zeros(numRuns,1);
r2TestAll    = zeros(numRuns,1);

tic;
numRowsTiles = ceil(sqrt(numRuns));
numColsTiles = ceil(numRuns/numRowsTiles);
figureHandle = figure('Position', [100, 100, numRowsTiles*800, numRowsTiles*600]);
t = tiledlayout(numRowsTiles, numColsTiles, 'Padding','compact','TileSpacing','compact');

for Run = 1:numRuns
    rng(Run*5);

    % shuffle rows
    nAll = size(X_mat,1);
    idx = randperm(nAll);
    Xs  = X_mat(idx,:);
    Ys  = Y_mat(idx,:);

    % splits
    trainRatio = 0.70; valRatio = 0.15; testRatio = 0.15;
    nTrain = floor(trainRatio*nAll);
    nVal   = floor(valRatio*nAll);
    idxTrain = 1:nTrain;
    idxVal   = (nTrain+1):(nTrain+nVal);
    idxTest  = (nTrain+nVal+1):nAll;

    X_train = Xs(idxTrain,:);  Y_train = Ys(idxTrain,:);
    X_val   = Xs(idxVal,:);    Y_val   = Ys(idxVal,:);
    X_test  = Xs(idxTest,:);   Y_test  = Ys(idxTest,:);

    % scale inputs (features only)
    mu    = mean(X_train,1);
    sigma = std(X_train,0,1);  sigma(sigma==0) = 1;
    X_train = (X_train - mu)./sigma;
    X_val   = (X_val   - mu)./sigma;
    X_test  = (X_test  - mu)./sigma;

    % -------------- Network --------------
    layers = [
        featureInputLayer(6,'Name','in')     % 6 features
        fullyConnectedLayer(128,'Name','fc1')
        reluLayer('Name','relu1')
        fullyConnectedLayer(64,'Name','fc2')
        reluLayer('Name','relu2')
        fullyConnectedLayer(5,'Name','out')  % 5 targets
        regressionLayer('Name','reg')
    ];

    miniB = 64;
    opts = trainingOptions('adam', ...
        'InitialLearnRate',0.0015, ...
        'MaxEpochs',10, ...
        'MiniBatchSize',miniB, ...,
        'Plots','training-progress',...
        'Shuffle','every-epoch', ...
        'ValidationData',{X_val,Y_val}, ...
        'ValidationFrequency',max(1,floor(size(X_train,1)/miniB)), ...
        'ValidationPatience',72, ...
        'L2Regularization',1e-6, ...
        'Verbose',true, ...
        'GradientThresholdMethod','l2norm', ...
        'GradientThreshold',1, ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropPeriod',100, ...
        'LearnRateDropFactor',0.84);

    % Train
    net = trainNetwork(X_train, Y_train, layers, opts);

    % Predict
    Y_train_pred = predict(net, X_train);
    Y_val_pred   = predict(net, X_val);
    Y_test_pred  = predict(net, X_test);

    % Metrics (aggregate over 5 outputs)
    Etr = Y_train_pred - Y_train;
    Eval= Y_val_pred   - Y_val;
    Ete = Y_test_pred  - Y_test;

    rmseTrainAll(Run) = aggRMSE(Etr);
    maeTrainAll(Run)  = aggMAE(Etr);
    r2TrainAll(Run)   = aggR2(Y_train_pred, Y_train);

    rmseValAll(Run) = aggRMSE(Eval);
    maeValAll(Run)  = aggMAE(Eval);
    r2ValAll(Run)   = aggR2(Y_val_pred, Y_val);

    rmseTestAll(Run) = aggRMSE(Ete);
    maeTestAll(Run)  = aggMAE(Ete);
    r2TestAll(Run)   = aggR2(Y_test_pred, Y_test);

    % Save model with metrics
    modelname = sprintf('Run_%d_Seed_%d_RMSE_%.4f_MAE_%.4f_R2_%.4f.mat', ...
        Run, Run*5, rmseTestAll(Run), maeTestAll(Run), r2TestAll(Run));
    save(fullfile(modelsDir, modelname), 'net','mu','sigma','featNames','targetCols');

    % Plot one target for sanity (m_H2_next, col 1)
    nexttile;
    scatter(Y_test(:,1), Y_test_pred(:,1), 10, 'filled'); hold on;
    limMax = 1.05*max([Y_test(:,1); Y_test_pred(:,1)]);
    plot([0 limMax],[0 limMax],'r--','LineWidth',1.2); hold off;
    xlabel('True m\_H2\_{next} [kg]'); ylabel('Predicted [kg]');
    title(sprintf('Run %d | R^2=%.3f', Run, r2TestAll(Run)));
    grid on;
end

% Save plots
figName = fullfile(plotsDir, 'AllRuns_PredVsTrue_mH2next.png');
exportgraphics(figureHandle, figName, 'Resolution', 300);

% ---------------- Metrics table ----------------
elapsedTime = toc;
fprintf('Elapsed Time        : %.2f seconds\n', elapsedTime);

metricsTable = table((1:numRuns)', ...
    rmseTrainAll, maeTrainAll, r2TrainAll, ...
    rmseValAll,   maeValAll,   r2ValAll, ...
    rmseTestAll,  maeTestAll,  r2TestAll, ...
    'VariableNames', {'Run', ...
    'RMSE_Train','MAE_Train','R2_Train', ...
    'RMSE_Val','MAE_Val','R2_Val', ...
    'RMSE_Test','MAE_Test','R2_Test'});

meanRow = {NaN, mean(rmseTrainAll), mean(maeTrainAll), mean(r2TrainAll), ...
                mean(rmseValAll),   mean(maeValAll),   mean(r2ValAll), ...
                mean(rmseTestAll),  mean(maeTestAll),  mean(r2TestAll)};

stdRow  = {NaN, std(rmseTrainAll), std(maeTrainAll), std(r2TrainAll), ...
                std(rmseValAll),   std(maeValAll),   std(r2ValAll), ...
                std(rmseTestAll),  std(maeTestAll),  std(r2TestAll)};

metricsTable = [metricsTable; ...
    cell2table(meanRow, 'VariableNames', metricsTable.Properties.VariableNames); ...
    cell2table(stdRow,  'VariableNames', metricsTable.Properties.VariableNames)];

metricsTable.Run = string(metricsTable.Run);
metricsTable.Run(end-1:end) = ["Mean"; "Std"];

disp(metricsTable);

csvName = fullfile(plotsDir, 'Metrics_AllRuns.csv');
writetable(metricsTable, csvName);
disp(['Metrics table saved as: ', csvName]);
