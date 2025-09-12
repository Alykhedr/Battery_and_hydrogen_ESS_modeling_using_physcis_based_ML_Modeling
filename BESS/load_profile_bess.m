function [time, P_el, P_load] = load_profile_bess()
    nHours = 365 * 24;
    time   = datetime(2023,1,1,0,0,0,'TimeZone','Europe/Berlin') + hours(0:nHours-1);

    % --- PV data import (retime to hourly) ---
    csvFile = 'Timeseries_52.522_13.360_SA3_1000kWp_crystSi_14_v45deg_2023_2023.csv';
    opts    = detectImportOptions(csvFile,'Delimiter',',','HeaderLines',10);
    T       = readtable(csvFile, opts);

    % Build a timetable with timezone-correct timestamps
    tUTC = datetime(T.time,'InputFormat','yyyyMMdd:HHmm','TimeZone','UTC');
    tBER = tUTC; tBER.TimeZone = 'Europe/Berlin';
    TT   = timetable(tBER, T.P, 'VariableNames', {'P'});

    % **Key step**: retime to hourly MEAN power
    TT_h = retime(TT, 'regular', 'mean', 'TimeStep', hours(1));

    % Align exactly 8760 hours
    if height(TT_h) < nHours
        error('PV series too short after hourly retime (%d < %d).', height(TT_h), nHours);
    end
    P_el = TT_h.P(1:nHours) / 1000;   % kW

    % --- Synthetic load (hourly) ---
    P_max_load = 500;                            % kW peak
    hoursDay   = hour(time) + minute(time)/60;
    doy        = day(time,'dayofyear');          % 1..365

    L1 = (sin(2*pi*(hoursDay-6)/24) + 1)/2;     % daily
    L2 = (sin(4*pi*(hoursDay-15)/24) + 1)/2;    % shoulders
    season_amp = 1 + 0.3*cos(2*pi*(doy-173)/365);

    P_load = P_max_load * season_amp .* (0.7*L1 + 0.3*L2) ...
           + 0.05*P_max_load*(rand(size(time))-0.5);
    P_load(P_load<0) = 0;
end
