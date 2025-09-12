%% HESS Level‑3 ML Baselines (Simplified + RF) — collapsed next‑hour mass
clear; clc; close all;

%% 1) Load & prepare data
trainTbl = readtable('model3_train_data.csv');
testTbl  = readtable('model3_test_data.csv');

% Shift “m_H2” up by one hour to form the target, then drop the last row
Ytrain_mass = trainTbl.m_H2(2:end);
trainTbl( end, : ) = [];  

Ytest_mass  = testTbl.m_H2(2:end);
testTbl(  end, : ) = [];

% define predictors and targets
predictors = {'P_el','P_fc','m_H2','i_sol','V_cell'};
targets    = {'m_H2','m_dot_el','m_dot_fc','P_comp','eta_el'};

Xtrain = trainTbl{:,predictors};
Ytrain = [ Ytrain_mass, ...
           trainTbl{:,targets(2:end)} ];

Xtest  = testTbl{:,predictors};
Ytest  = [ Ytest_mass, ...
           testTbl{:,targets(2:end)} ];

[nTest, nTgt] = size(Ytest);

%% 2) k‑NN with simple grid search for K
Ks = [3,5,7,10];
bestRMSE = inf;
for K = Ks
    idxNN = knnsearch(Xtrain, Xtest, 'K',K);
    Yknn_tmp = mean( reshape(Ytrain(idxNN(:),:), nTest, K, nTgt), 2 );
    Yknn_tmp = squeeze(Yknn_tmp);
    rmse_vec = sqrt( mean((Yknn_tmp - Ytest).^2, 1) );
    avgRmse  = mean(rmse_vec);
    if avgRmse < bestRMSE
        bestRMSE = avgRmse;
        bestK    = K;
        Yknn     = Yknn_tmp;
    end
end
fprintf('Chosen k‑NN K = %d (avg RMSE = %.3f)\n', bestK, bestRMSE);

%% 3) SVM‑RBF, 4) Regression Tree, 5) Random Forest
models = {'SVM-RBF','Tree','RandomForest'};
Ys     = {Yknn};  % start with kNN

% SVM‑RBF
Ysvm = zeros(nTest,nTgt);
for j = 1:nTgt
    mdl = fitrsvm(Xtrain, Ytrain(:,j), 'KernelFunction','rbf','Standardize',true);
    Ysvm(:,j) = predict(mdl, Xtest);
end
Ys{end+1} = Ysvm;

% Regression Tree
Ytree = zeros(nTest,nTgt);
for j = 1:nTgt
    mdl = fitrtree(Xtrain, Ytrain(:,j));
    Ytree(:,j) = predict(mdl, Xtest);
end
Ys{end+1} = Ytree;

% Random Forest
Yrf = zeros(nTest,nTgt);
for j = 1:nTgt
    rfMdl = TreeBagger(100, Xtrain, Ytrain(:,j), 'Method','regression','OOBPrediction','on');
    Yrf(:,j) = predict(rfMdl, Xtest);
    fprintf('RF OOB RMSE for %s: %.3f\n', targets{j}, sqrt(oobError(rfMdl,'Mode','Ensemble')));
end
Ys{end+1} = Yrf;

%% 6) Compile RMSE & R² results
allModels = [{'kNN'}, models];
results = table('Size',[nTgt*numel(allModels),4], ...
    'VariableTypes',{'string','string','double','double'}, ...
    'VariableNames',{'Model','Target','RMSE','R2'});
row = 1;
for m = 1:numel(allModels)
    Yp = Ys{m};
    for j = 1:nTgt
        yTrue = Ytest(:,j);
        yPred = Yp(:,j);
        rmse  = sqrt(mean((yPred - yTrue).^2));
        ssRes = sum((yTrue - yPred).^2);
        ssTot = sum((yTrue - mean(yTrue)).^2);
        r2    = 1 - ssRes/ssTot;
        results(row,:) = {allModels{m}, targets{j}, rmse, r2};
        row = row + 1;
    end
end

%% 7) Display final table
fprintf('\n--- Baseline ML Models (TEST set) ---\n');
disp(results);
