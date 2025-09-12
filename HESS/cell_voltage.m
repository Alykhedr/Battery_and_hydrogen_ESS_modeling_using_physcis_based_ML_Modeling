% File: hess_helpers.m
% Common helper functions for HESS models (3 & 4)

function V = cell_voltage(i, p)
% CELL_VOLTAGE  Calculates a single-cell voltage [V]
%   i [A/cm^2]   : current density
%   p            : parameter struct from params.m

  F    = p.F;               % Faraday constant [C/mol]
  T    = p.el.T;         % Electrolyzer temperature [K]
  PH2  = p.el.P_H2;         % H2 partial pressure [bar]
  PO2  = p.el.P_O2;         % O2 partial pressure [bar]
  beta = p.el.beta_H2O;     % water activity (≈1)

  % Reversible voltage term
  Vrev = p.el.V_rev + p.R*T/(2*F)*log(PH2*sqrt(PO2)/beta);
  % Ohmic overpotential
  Vohm = p.el.Ohmic * i;
  % Activation overpotential
  Vact = p.R*T/(p.el.sigma_p*F)*asinh(i/(2*p.el.i_p)) + ...
         p.R*T/(p.el.sigma_n*F)*asinh(i/(2*p.el.i_n));

  V = Vrev + Vohm + Vact;
end

