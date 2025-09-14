
clear; clc;

%% ---------- SETTINGS ----------
inputFile     = 'DischargeCapacity.xlsx';
SCALE_METHOD  = 'standard';   % 'standard' | 'minmax' | 'none'

cleanCSV      = 'calendar_clean.csv';
normCSV       = 'calendar_clean_norm.csv';

%% ---------- STEP 1: STRUCTURE (wide -> long) ----------
T = readtable(inputFile, 'VariableNamingRule','preserve');

% -- locate time column --
vnames = lower(string(T.Properties.VariableNames));
cand   = ["t_h","time","t","hours"];
idxT   = find(ismember(vnames, cand), 1);
if isempty(idxT), idxT = 1; end

% -- time to hours --
t_h = T{:, idxT};
if isdatetime(t_h)
    t_h = hours(t_h - t_h(1));
else
    t_h = double(t_h);
end
nT = height(T);

% -- condition columns & total rows --
condCols = setdiff(1:width(T), idxT);
nC = numel(condCols);
N  = nT * nC;

% -- PREALLOCATE (no growth in loop) --
cond_id_all   = strings(N,1);
T_degC_all    = zeros(N,1);
T_K_all       = zeros(N,1);
SOC_pct_all   = zeros(N,1);
SOC_store_all = zeros(N,1);
t_h_all       = zeros(N,1);
Q_Ah_all      = zeros(N,1);
source_all    = strings(N,1);

% -- fill by blocks --
for k = 1:nC
    c   = condCols(k);
    hdr = string(T.Properties.VariableNames{c});

    % parse T (°C) and SOC (%)
    tokT = regexp(hdr, '(-?\d+\.?\d*)\s*°?\s*C', 'tokens', 'once');
    assert(~isempty(tokT), "Cannot parse temperature: " + hdr);
    T_degC = str2double(tokT{1});  T_K = T_degC + 273.15;

    tokS = regexp(hdr, '(\d+\.?\d*)\s*%?\s*SOC', 'tokens', 'once');
    assert(~isempty(tokS), "Cannot parse SOC: " + hdr);
    SOC_pct   = str2double(tokS{1});
    SOC_store = SOC_pct/100;

    cond_id = sprintf('TP_%sC_%sSOC', erase(sprintf('%.0f',T_degC),'+'), sprintf('%.0f',SOC_pct));

    q = double(T{:, c});  % capacity column (Q_Ah)

    % target rows for this condition
    rows = (k-1)*nT + (1:nT);

    % assign
    cond_id_all(rows)   = cond_id;
    T_degC_all(rows)    = T_degC;
    T_K_all(rows)       = T_K;
    SOC_pct_all(rows)   = SOC_pct;
    SOC_store_all(rows) = SOC_store;
    t_h_all(rows)       = t_h;
    Q_Ah_all(rows)      = q;
    source_all(rows)    = hdr;
end

out = table( ...
    cond_id_all, T_degC_all, T_K_all, SOC_store_all, SOC_pct_all, ...
    t_h_all, Q_Ah_all, source_all, ...
    'VariableNames', {'cond_id','T_degC','T_K','SOC_store','SOC_pct','t_h','Q_Ah','source_col'});

writetable(out, cleanCSV);
disp("Wrote: " + cleanCSV + "  (#rows=" + height(out) + ")");

%% ---------- STEP 2: MINIMAL CLEAN + NORMALIZE ----------

% Assumptions: no outliers, no negatives; we keep Q_Ah raw.
T2 = sortrows(out, {'cond_id','t_h'});

% -- SOH_init: vs first valid Q_Ah per condition --
[grp, ids] = findgroups(T2.cond_id);
Q0 = splitapply(@(q) q(find(~isnan(q),1,'first')), T2.Q_Ah, grp);
mQ0 = containers.Map(ids, num2cell(Q0));
T2.SOH_init = arrayfun(@(i) T2.Q_Ah(i) / mQ0(T2.cond_id{i}), (1:height(T2))');


% -- EFC placeholder (calendar-only) --
if ~ismember('EFC', string(T2.Properties.VariableNames))
    T2.EFC = zeros(height(T2),1);
end

% -- Normalize inputs (T_degC, SOC_store, t_h) --
Xvars = {'T_degC','SOC_store','t_h'};
X = T2{:, Xvars};

switch lower(SCALE_METHOD)
    case 'standard'
        mu = mean(X,1,'omitnan'); sg = std(X,0,1,'omitnan'); sg(sg==0)=1;
        Xn = (X - mu) ./ sg;
    case 'minmax'
        mn = min(X,[],1); mx = max(X,[],1); span = mx - mn; span(span==0)=1;
        Xn = (X - mn) ./ span;
    case 'none'
        Xn = X;
    otherwise
        error('SCALE_METHOD must be ''standard'', ''minmax'', or ''none''.');
end

for j = 1:numel(Xvars)
    T2.([Xvars{j} '_n']) = Xn(:, j);
end

% -- Final column order --
T2 = T2(:, {'cond_id','T_degC','T_K','SOC_store','SOC_pct','t_h','Q_Ah', ...
            'SOH_init','EFC','source_col', ...
            'T_degC_n','SOC_store_n','t_h_n'});

writetable(T2, normCSV);
disp("Wrote: " + normCSV + "  (#rows=" + height(T2) + ", conds=" + numel(unique(T2.cond_id)) + ")");

%% ---------- STEP 3: FEATURE ENGINEERING & DATASET SPLITTING (FFNN) ----------
% Goal:
%   - Inputs (X): T_degC, SOC_store, t_h (scaled using TRAIN-only stats)
%   - Target (y): SOH = SOH_init (per-condition BOL normalization)
%   - Stratified split by cond_id so every condition appears in all splits (if possible)

rng(42);  % reproducible shuffle

% ----- 3A. Choose features/target from the unscaled columns -----
featNames = {'T_degC','SOC_store','t_h'};   % raw features
targetCol = 'SOH_init';                     % target (0..1)

% We'll compute scalers using TRAIN only (ignoring the earlier dataset-level scaling).
X_all = T2{:, featNames};
y_all = T2{:, targetCol};

% ----- 3B. Stratified split by condition -----
u = unique(T2.cond_id, 'stable');

trainCells = cell(numel(u),1);
valCells   = cell(numel(u),1);
testCells  = cell(numel(u),1);

for k = 1:numel(u)
    idx = find(strcmp(T2.cond_id, u{k}));
    idx = idx(randperm(numel(idx)));  % shuffle within condition
    n   = numel(idx);

    if n >= 3
        nTr = max(1, round(0.70*n));
        nVa = max(1, round(0.15*n));
        nTe = n - nTr - nVa;
        if nTe < 1
            nTe = 1; 
            nVa = max(1, n - nTr - nTe);
        end
        takeTr = idx(1:nTr);
        takeVa = idx(nTr+1:nTr+nVa);
        takeTe = idx(nTr+nVa+1:end);
    elseif n == 2
        takeTr = idx(1);
        takeVa = []; 
        takeTe = idx(2);
    else % n == 1
        takeTr = idx(1);
        takeVa = []; 
        takeTe = [];
    end

    trainCells{k} = takeTr(:);
    valCells{k}   = takeVa(:);
    testCells{k}  = takeTe(:);
end

% Concatenate once (no size changes in loop)
idxTrain = vertcat(trainCells{:});
idxVal   = vertcat(valCells{:});
idxTest  = vertcat(testCells{:});

% Optional: randomize order within each split
if ~isempty(idxTrain), idxTrain = idxTrain(randperm(numel(idxTrain))); end
if ~isempty(idxVal),   idxVal   = idxVal(randperm(numel(idxVal)));   end
if ~isempty(idxTest),  idxTest  = idxTest(randperm(numel(idxTest)));  end

% ----- 3C. Build train/val/test tables -----
TrainTbl = T2(idxTrain, :);
ValTbl   = T2(idxVal,   :);
TestTbl  = T2(idxTest,  :);

% Sanity check (counts)
assert(height(TrainTbl) + height(ValTbl) + height(TestTbl) == height(T2), 'Split counts mismatch.');

% Ensure each condition appears in train (by construction it should)
% (Val/Test coverage depends on per-condition sample count.)

% ----- Compute scalers on TRAIN ONLY -----
Xtr = TrainTbl{:, featNames};

switch lower(SCALE_METHOD)
    case 'standard'
        mu = mean(Xtr, 1, 'omitnan');
        sg = std(Xtr, 0, 1, 'omitnan'); sg(sg==0) = 1;
        scaleFun  = @(X) (X - mu) ./ sg;
        scalers = struct('method','standard','feature_names',{featNames},'mu',mu,'sigma',sg);

    case 'minmax'
        mn = min(Xtr, [], 1);
        mx = max(Xtr, [], 1);
        span = mx - mn; span(span==0) = 1;
        scaleFun  = @(X) (X - mn) ./ span;
        scalers = struct('method','minmax','feature_names',{featNames},'min',mn,'max',mx);

    case 'none'
        scaleFun  = @(X) X;
        scalers = struct('method','none','feature_names',{featNames});

    otherwise
        error('SCALE_METHOD must be ''standard'', ''minmax'', or ''none''.');
end

% Apply scaler to all splits
X_train = scaleFun(TrainTbl{:, featNames});
X_val   = scaleFun(ValTbl{:,   featNames});
X_test  = scaleFun(TestTbl{:,  featNames});

y_train = TrainTbl{:, targetCol};
y_val   = ValTbl{:,   targetCol};
y_test  = TestTbl{:,  targetCol};

% (Optional) attach scaled columns to tables for inspection
for j = 1:numel(featNames)
    TrainTbl.([featNames{j} '_trn']) = X_train(:, j);
    ValTbl.([featNames{j} '_trn'])   = X_val(:, j);
    TestTbl.([featNames{j} '_trn'])  = X_test(:, j);
end

% ----- Quick verification of splits -----
fprintf('Split sizes: train=%d, val=%d, test=%d (total=%d)\n', ...
    numel(y_train), numel(y_val), numel(y_test), height(T2));

% Distribution check (means & stds) — TRAIN vs VAL/TEST (feature-wise)
if ~strcmpi(SCALE_METHOD,'none')
    mTrain = mean(X_train,1); sTrain = std(X_train,0,1);
    mVal   = mean(X_val,1);   sVal   = std(X_val,0,1);
    mTest  = mean(X_test,1);  sTest  = std(X_test,0,1);
    disp(table(featNames', mTrain', sTrain', mVal', sVal', mTest', sTest', ...
        'VariableNames', {'Feature','TrainMean','TrainStd','ValMean','ValStd','TestMean','TestStd'}));
end

% Optional: ensure some conditions appear across splits (info only)
% (Small n per condition may prevent full coverage.)
% disp(groupsummary(TrainTbl,'cond_id'));
% disp(groupsummary(ValTbl,'cond_id'));
% disp(groupsummary(TestTbl,'cond_id'));

% ----- Save artifacts -----
save('calendar_splits.mat', ...
     'featNames','targetCol','scalers', ...
     'X_train','y_train','X_val','y_val','X_test','y_test', ...
     'TrainTbl','ValTbl','TestTbl');

writetable(TrainTbl, 'calendar_train.csv');
if ~isempty(ValTbl),  writetable(ValTbl,  'calendar_val.csv');  end
if ~isempty(TestTbl), writetable(TestTbl, 'calendar_test.csv'); end

disp('Saved: calendar_splits.mat + CSVs (train/val/test). Ready for Step 4 (FFNN).');

%% ---------- STEP 4: FFNN DESIGN & TRAINING (Calendar Aging) ----------
% Use 2-D arrays: X_* are [N x 3], y_* are [N x 1]
XTrain = single(X_train);   YTrain = single(y_train(:));
XVal   = single(X_val);     YVal   = single(y_val(:));
XTest  = single(X_test);    YTest  = single(y_test(:));

% 4A) Network (3 -> 32 -> 16 -> 1)
layers = [
    featureInputLayer(3,'Name','in')
    fullyConnectedLayer(32,'Name','fc1')
    reluLayer('Name','relu1')
    dropoutLayer(0.1,'Name','drop1')
    fullyConnectedLayer(16,'Name','fc2')
    reluLayer('Name','relu2')
    fullyConnectedLayer(1,'Name','out')
    regressionLayer('Name','reg')
];

% 4B) Training options (val data as 2-D arrays)
miniB = 32;
opts = trainingOptions('adam', ...
    'InitialLearnRate',1e-2, ...
    'MaxEpochs',200, ...
    'MiniBatchSize',miniB, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XVal, YVal}, ...
    'ValidationFrequency',max(1, floor(size(XTrain,1)/miniB)), ...
    'ValidationPatience',10, ...
    'L2Regularization',1e-4, ...
    'Verbose',false, ...
    'Plots','training-progress','GradientThresholdMethod','l2norm','GradientThreshold',1, ...
'LearnRateSchedule','piecewise','LearnRateDropPeriod',10,'LearnRateDropFactor',0.5);

rng(1);
net = trainNetwork(XTrain, YTrain, layers, opts);

% 4C) Evaluate
y_pred_train = predict(net, XTrain);  y_pred_train = gather(y_pred_train(:));
y_pred_val   = predict(net, XVal);    y_pred_val   = gather(y_pred_val(:));
y_pred_test  = predict(net, XTest);   y_pred_test  = gather(y_pred_test(:));

% Metrics
mse  = @(a,b) mean((a-b).^2);
rmse = @(a,b) sqrt(mse(a,b));
mae  = @(a,b) mean(abs(a-b));
r2   = @(a,b) 1 - sum((a-b).^2)/sum((b-mean(b)).^2);

RMSE = [ rmse(y_pred_train,y_train);
         rmse(y_pred_val,  y_val);
         rmse(y_pred_test, y_test) ];

MAE  = [ mae(y_pred_train,y_train);
         mae(y_pred_val,  y_val);
         mae(y_pred_test, y_test) ];

R2   = [ r2(y_pred_train,y_train);
         r2(y_pred_val,  y_val);
         r2(y_pred_test, y_test) ];

metrics = table(RMSE, MAE, R2, 'RowNames', {'Train','Val','Test'});
disp(metrics)

% 4C) Plot 1: Pred vs True (Test)
figure('Name','Pred vs True (Test)'); grid on; hold on;
scatter(y_test, y_pred_test, 20, 'filled');
plot([min(y_test) max(y_test)], [min(y_test) max(y_test)], 'k--', 'LineWidth',1);
xlabel('True SOH'); ylabel('Predicted SOH'); title('FFNN — Test Set'); axis tight;
