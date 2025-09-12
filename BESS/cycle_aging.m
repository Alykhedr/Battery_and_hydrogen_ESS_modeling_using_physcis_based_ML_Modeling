function delta = cycle_aging(Cch, Cdis, mSOC, p, fT, fD)
% MONTES_DELTA  Percent-points per EFC^a for Montes LFP cycle law.

delta = p.kcyc * fT * fD * exp(p.kCch*Cch) * exp(p.kCdch*Cdis) ...
      * (1 + p.kmSOC*mSOC*((1 - mSOC)/(2*p.mSOCref)));
end
