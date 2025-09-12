function V = ocv_from_soc(SOC)
%OCV_FROM_SOC Returns open-circuit voltage [V] for given SOC
%  Validated Samsung 21700 (INR21700) OCV curve (16 points)
  soc_table = [0.0000 0.0236 0.0473 0.0709 0.0945 0.1238 0.1530 0.2417 ...
               0.3303 0.4644 0.5985 0.7391 0.8798 0.9199 0.9599 1.0000];
  ocv_table = [2.6929 3.1683 3.3177 3.3668 3.3923 3.4225 3.4561 3.5478 ...
               3.6094 3.7059 3.8368 3.9740 4.0759 4.1018 4.1315 4.1710];
  % Piecewise cubic interpolation for smooth curve
  V = interp1(soc_table, ocv_table, SOC, 'pchip');
end