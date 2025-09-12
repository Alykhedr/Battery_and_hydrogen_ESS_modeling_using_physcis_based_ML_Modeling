function [dQcal_pct, fTcal, fSOC] = calendar_aging(SOC, dt_h, tcal_h, p, fTcal_const)
% CALENDAR_INCREMENT  Percent-points capacity loss for calendar over dt_h.
% Uses Naumann/Schimpe square-root time law with Arrhenius + anode-OCV SOC dependence.

% Temperature factor
if nargin >= 5 && ~isempty(fTcal_const)
    fTcal = fTcal_const;
else
    fTcal = exp( -p.deg.Ea_cal_Jmol/p.R_gas * (1/(p.T_degC+273.15) - 1/p.deg.Tref_K) );
end

% SOC dependence (anode OCV → Tafel-like multiplier)
if p.modes.use_soc_dependence && isfield(p.deg,'Ua_ref_V') && ~isnan(p.deg.Ua_ref_V)
    Ua   = Ua_graphite(max(0,min(1,SOC)), p.deg.OCVmap);
    fSOC = ( exp((p.deg.alpha*p.F/p.R_gas) * ((p.deg.Ua_ref_V - Ua)/p.deg.Tref_K)) + p.deg.k0 ) ...
           / (1 + p.deg.k0);
else
    fSOC = 1.0;
end

% Square-root in time (convert to percent points)
dQcal_pct = 100 * p.deg.kcal_ref * fTcal * fSOC * ( sqrt(tcal_h + dt_h) - sqrt(tcal_h) );
end
