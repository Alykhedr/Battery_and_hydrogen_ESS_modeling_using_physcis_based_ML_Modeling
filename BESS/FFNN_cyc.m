% cycle_step1_build_long_min.m
% STEP 1 ONLY — Parse & organize cycle Capacity–FEC data (wide -> long), minimal version.
% Includes any *.mat whose name contains "Capacity_CC_CV_FEC".
% Excludes *_Time* (we're building the FEC table).
% Outputs:
%   - cycle_clean.csv
%   - cycle_clean.mat  (T_cyc)

clear; clc;

outCSV = 'cycle_clean.csv';
outMAT = 'cycle_clean.mat';

%% -------- Collect files (simple filter; R files already removed) --------
D = dir('*Capacity_CC_CV_FEC*.mat');
% drop time-based capacity files (keep only FEC)
mask = ~contains({D.name}, '_Time', 'IgnoreCase', true);
D = D(mask);
files = fullfile({D.folder}, {D.name});
assert(~isempty(files), 'No Capacity_CC_CV_FEC .mat files found.');

%% -------- Pass 1: count rows (for preallocation) --------
totalRows = 0;
for f = 1:numel(files)
    S = load(files{f});
    if ~isfield(S,'X_Axis_Data_Mat') || ~isfield(S,'Y_Axis_Data_Mat'), continue; end
    X = S.X_Axis_Data_Mat; Y = S.Y_Axis_Data_Mat;
    if ~isequal(size(X), size(Y)), continue; end
    [~,N] = size(X);
    for j = 1:N
        m = ~isnan(X(:,j)) & ~isnan(Y(:,j));
        c = nnz(m);
        if c > 0, totalRows = totalRows + c; end
    end
end
assert(totalRows > 0, 'No finite (FEC, SOH) pairs found.');

%% -------- Preallocate essentials --------
proto_id  = strings(totalRows,1);
FEC       = zeros(totalRows,1);
SOH       = zeros(totalRows,1);
T_degC    = NaN(totalRows,1);
mean_SOC  = NaN(totalRows,1);   % (%)
DOD_pct   = NaN(totalRows,1);   % (%)
Cch_rate  = NaN(totalRows,1);
Cdch_rate = NaN(totalRows,1);

% helpers
parseLegend = @(s) local_parse_legend(s);
cleanId     = @(s) regexprep(regexprep(regexprep(regexprep(string(s),'\s+','_'),'°','deg'),'%','pct'),'[^A-Za-z0-9_\-\.]','');

%% -------- Pass 2: fill arrays --------
p = 1;
for f = 1:numel(files)
    S = load(files{f});
    if ~isfield(S,'X_Axis_Data_Mat') || ~isfield(S,'Y_Axis_Data_Mat'), continue; end
    X = S.X_Axis_Data_Mat; Y = S.Y_Axis_Data_Mat;
    [~,N] = size(X);

    if isfield(S,'Legend_Vec') && numel(S.Legend_Vec)==N
        L = S.Legend_Vec;
    else
        L = repmat({''},1,N);
    end

    for j = 1:N
        m = ~isnan(X(:,j)) & ~isnan(Y(:,j));
        if ~any(m), continue; end

        idx  = find(m);
        n    = numel(idx);
        rows = p:(p+n-1);

        % data
        FEC(rows) = X(idx,j);
        SOH(rows) = Y(idx,j);

        % meta (parsed from legend)
        meta = parseLegend(L{j});
        T_degC(rows)    = meta.T_degC;
        mean_SOC(rows)  = meta.mean_SOC;
        DOD_pct(rows)   = meta.DOD_pct;
        Cch_rate(rows)  = meta.Cch_rate;
        Cdch_rate(rows) = meta.Cdch_rate;

        % proto id
        pid = cleanId(L{j});
        if strlength(pid) == 0
            pid = cleanId(sprintf('Proto_T%g_SOC%g_DOD%g_C+%g_C-%g', ...
                  meta.T_degC, meta.mean_SOC, meta.DOD_pct, meta.Cch_rate, meta.Cdch_rate));
        end
        proto_id(rows) = pid;

        p = p + n;
    end
end

% trim if over-allocated
if p <= totalRows
    keep = 1:(p-1);
    proto_id  = proto_id(keep);
    FEC       = FEC(keep);
    SOH       = SOH(keep);
    T_degC    = T_degC(keep);
    mean_SOC  = mean_SOC(keep);
    DOD_pct   = DOD_pct(keep);
    Cch_rate  = Cch_rate(keep);
    Cdch_rate = Cdch_rate(keep);
end

%% -------- Build table, sort, save --------
T_cyc = table(proto_id, FEC, SOH, T_degC, mean_SOC, DOD_pct, Cch_rate, Cdch_rate);
T_cyc = sortrows(T_cyc, {'proto_id','FEC'});

writetable(T_cyc, outCSV);
save(outMAT, 'T_cyc');

%% ================= local parser =================
function meta = local_parse_legend(lbl)
% Parse legend strings like:
%  'Testpoint Cyclization_40°C_50%SOC_100%DOD_1C_1C_CC+CV'
%  'Testpoint LoadSpectrumPVBattery_40°C_51.4%SOC'
% Missing tokens -> NaN (numeric), "" (string ignored here)

    s = string(lbl);

    % Temperature (°C)
    T = NaN; m = regexp(s, '(-?\d+(?:\.\d+)?)\s*°C', 'tokens', 'once');
    if ~isempty(m), T = str2double(m{1}); end

    % Mean SOC (%)
    mSOC = NaN; m = regexp(s, '(\d+(?:\.\d+)?)\s*%SOC', 'tokens', 'once');
    if ~isempty(m), mSOC = str2double(m{1}); end

    % DOD (%)
    DOD = NaN; m = regexp(s, '(\d+(?:\.\d+)?)\s*%DOD', 'tokens', 'once');
    if ~isempty(m), DOD = str2double(m{1}); end

    % C-rates "xC_yC"
    Cch = NaN; Cdch = NaN;
    m = regexp(s, '(\d+(?:\.\d+)?)C_(\d+(?:\.\d+)?)C', 'tokens', 'once');
    if ~isempty(m)
        Cch  = str2double(m{1});
        Cdch = str2double(m{2});
    end

    meta = struct('T_degC',T,'mean_SOC',mSOC,'DOD_pct',DOD,'Cch_rate',Cch,'Cdch_rate',Cdch);
end
