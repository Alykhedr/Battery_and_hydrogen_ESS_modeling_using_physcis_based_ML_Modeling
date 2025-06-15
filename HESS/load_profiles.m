function [time, P_el, P_fc] = load_profiles(Timeseries)
  % expects a table with columns: time, P_el, P_fc
  T = readtable(Timeseries);
  time = T.time;
  P_el = T.P;
  % if your file has G(i) or H_sun etc, pick whichever column(s) map to P_el/P_fc
  P_fc = zeros(size(P_el));  % or load a second column if you have FC profiles
end