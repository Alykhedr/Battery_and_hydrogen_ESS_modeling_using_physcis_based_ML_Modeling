function p = params(time)
  % — Constants —
  p.LHV_H2      = 33.33;   % kWh/kg
  p.eta_el      = 0.75;
  p.eta_fc      = 0.60;
  p.max_mass_H2 = 50;      % kg

  % — Pre-allocate state arrays —
  N = numel(time);
  p.init.m_dot_el = zeros(1, N);
  p.init.m_dot_fc = zeros(1, N);
  p.init.mass_H2  = zeros(1, N);

  % — Initial tank level —
  p.init.mass_H2(1) = 0.5 * p.max_mass_H2;
end
