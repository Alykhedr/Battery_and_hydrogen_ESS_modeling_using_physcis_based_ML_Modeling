function Vcell = fc_voltage(j,p)
  % j: current density [A/cm^2]
  R = p.R; F = p.F; T = p.FC.T;
  PH2 = p.FC.p_H2;  PO2 = p.FC.p_O2;

  % 1) Nernst
  ENer = 1.229 ...
        - 8.5e-4*(T - 298.15) ...
        + 4.308e-5*T*(log(PH2) + 0.5*log(PO2));

  % 2) activation
  z = p.FC.zeta;
  Eact = z(1) + z(2)*T + z(3)*T*log(j) + z(4)*T*log(j);

  % 3) ohmic (membrane + contact)
  rhoM = p.FC.rho0*(1 + p.FC.rhoa*j + p.FC.rhob*(j.^2.5));
  RM   = rhoM * (p.el.L*1e-2)/p.FC.A;    % cm→m convert path length
  Rohm = j*(RM + p.FC.Rc);
  Eohm = Rohm;

  % 4) concentration
  Econ = (R*T/(2*F))*log(1 - j/p.FC.j_lim);

  Vcell = ENer - Eact - Eohm - Econ;
end
