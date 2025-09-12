function [time, P_el, P_fc] = load_profiles_2()
%LOAD_PROFILES_2  168-h PV & FC profiles for August week
%   Real PV from PVGIS CSV + synthetic FC demand proxy based on ENTSO-E (2022) data
%   PV clear-sky August Berlin 1 MWp SARAH3 citePVGIS2023CSV
%   FC demand proxy: German grid-load shape from ENTSO-E (2022) citeENTSOE2022Load

  %% 1) Time axis
  time = (0:167)';  % hours since start of week

  %% 2) PV input: first 168 hours of August from PVGIS CSV
  csvFile = 'Timeseries_52.522_13.360_SA3_1000kWp_crystSi_14_v45deg_2023_2023.csv';
  opts = detectImportOptions(csvFile, 'Delimiter', ',', 'HeaderLines', 10);
  T = readtable(csvFile, opts);
  % parse and filter for August
  T.datetime = datetime(T.time, 'InputFormat', 'yyyyMMdd:HHmm', 'TimeZone','UTC');
  T.datetime.TimeZone = 'Europe/Berlin';
  T_aug = T(month(T.datetime)==8, :);
  assert(height(T_aug)>=168, 'Not enough August data');
  T_week = T_aug(1:168, :);
  % convert P from W→kW
  P_el = T_week.P / 1000;
  assert(numel(P_el)==168 && all(P_el>=0), 'Invalid PV profile');

      %% 3) FC demand: fuel-cell refuelling station profile (synthetic, event-based)
  base = 10;                 % kW standby
  P_fc = base * ones(168,1);
  rng(99);                   % reproducible events
  events_per_day = 3;
  for d = 0:6
    hours = randperm(13, events_per_day) + 6 + 24*d;  % choose hours 7–19
    for h = hours
      amp = 150 + 150*rand;  % [150–300] kW
      idx = h + 1;           % MATLAB index
      if idx <= 168
        P_fc(idx) = amp;
      end
    end
  end
  P_fc = P_fc + 5*randn(168,1); % add noise
  P_fc(P_fc < 0) = 0;
  P_fc(P_fc > 200) = 200;        % cap at FC rating
  assert(numel(P_fc)==168 && all(P_fc>=0), 'Invalid FC profile');
end
