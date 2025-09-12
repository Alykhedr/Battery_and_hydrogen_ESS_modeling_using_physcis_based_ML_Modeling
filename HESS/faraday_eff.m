function etaF = faraday_eff(i_el, p)
% FARADAY_EFFICIENCY  Computes Faraday efficiency [-]
%   i_el [A/cm^2] : current density
%   p             : parameter struct from params.m
  F    = p.F;         % Faraday constant [C/mol]
  L    = p.el.L;      % Membrane thickness [cm]
  % Permeability coefficients [mol/(cm·s·bar)]
  eH2f = p.el.epsilon_H2f;
  eO2f = p.el.epsilon_O2f;
  eH2p = p.el.epsilon_H2p;
  PH2  = p.el.P_H2;   % H2 pressure [bar]
  PO2  = p.el.P_O2;   % O2 pressure [bar]

  % Combined crossover term (A/cm^2)
  j_cross = (2*F*eH2f*PH2 + 4*F*eO2f*PO2 + 2*F*eH2p*PH2*PO2) / L;
  eta_raw = 1 - j_cross ./ max(i_el, 1e-4);  
  % Faraday efficiency (unitless)
  etaF   = min(1, max(0, eta_raw));  % clamp into [0,1]
end
