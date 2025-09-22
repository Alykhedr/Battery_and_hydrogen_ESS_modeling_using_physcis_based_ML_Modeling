%% build_PV_FC_profiles_from_DE.m
% Input:  CSV with columns:
%   utc_timestamp, DE_load_actual_entsoe_transparency, DE_solar_capacity, DE_solar_generation_actual
% Output: CSV with: time, P_pv (kW), P_load (kW)

infile  = 'DE_2018Apr_2020Jan_60min.csv';   % <-- your generated file
outfile = 'profiles_PV1MW_FCscaled_60min.csv';

% ---- user knobs ----
Ppv_rated_MW   = 3.0;   % target PV plant size (MW)
Pfc_rated_kW   = 400;   % your FC rating (kW)
u_target       = 0.30;  % desired average FC utilization (0..1)
clip_to_rating = true;  % cap P_load at Pfc_rated_kW (recommended)

% ---- read (keep timestamp as text; we’ll parse ourselves) ----
opts = detectImportOptions(infile,'NumHeaderLines',0);
opts = setvaropts(opts, 'utc_timestamp', 'Type', 'char');
T = readtable(infile, opts);

% sanity names (handle possible slight header diffs)
must = ["utc_timestamp","DE_load_actual_entsoe_transparency","DE_solar_capacity","DE_solar_generation_actual"];
assert(all(ismember(must, string(T.Properties.VariableNames))), 'Input file missing required columns.');

% ---- parse time (UTC) ----
ts = string(T.utc_timestamp);
t  = datetime(ts, 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'", 'TimeZone','UTC');

% ---- clean numeric series (fill short gaps) ----
load_DE = T.DE_load_actual_entsoe_transparency;      % MW
cap_pv  = T.DE_solar_capacity;                        % MW
gen_pv  = T.DE_solar_generation_actual;               % MW

% Fill ≤6h gaps using time-aware interpolation
for v = 1:3
    x = {load_DE, cap_pv, gen_pv};
    xi = fillmissing(x{v}, 'linear', 'SamplePoints', t, 'MaxGap', hours(6));
    switch v
        case 1, load_DE = xi;
        case 2, cap_pv  = xi;
        case 3, gen_pv  = xi;
    end
end

% ---- PV to 1 MW via capacity factor ----
% CF = generation / capacity (both MW). If capacity <= 0 → CF = 0.
CF = zeros(size(gen_pv));
mask_cap = cap_pv > 0;
CF(mask_cap) = gen_pv(mask_cap) ./ cap_pv(mask_cap);
% numeric hygiene: clip tiny negatives / spikes
CF = max(CF, 0);
CF = min(CF, 1.2);   % allow slight >1 due to reporting noise, but tame it

P_pv_kW = CF * Ppv_rated_MW * 1000;     % kW

% ---- Load scaled to FC utilization target ----
% We want mean(P_load) ≈ u_target * Pfc_rated_kW
mean_nat_MW = mean(load_DE, 'omitnan');     % MW (national)
target_mean_kW = u_target * Pfc_rated_kW;   % kW
scale = (target_mean_kW/1000) / mean_nat_MW; % (kW→MW) / MW = scalar

P_load_kW = load_DE * scale * 1000;   % kW

if clip_to_rating
    P_load_kW = min(P_load_kW, Pfc_rated_kW);
    P_load_kW = max(P_load_kW, 0);
end

% ---- emit tidy table ----
Tout = table( ...
    string(t, "yyyy-MM-dd'T'HH:mm:ss'Z'"), ...
    round(P_pv_kW, 3), ...
    round(P_load_kW, 3), ...
    'VariableNames', {'time','P_pv','P_load'} ...
);

writetable(Tout, outfile);
fprintf('Wrote %s  | rows=%d  | avg P_pv=%.1f kW  | avg P_load=%.1f kW (target %.1f kW)\n', ...
        outfile, height(Tout), mean(Tout.P_pv), mean(Tout.P_load), u_target*Pfc_rated_kW);
