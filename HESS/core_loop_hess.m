
function out = core_loop_hess(time, P_el, P_fc, p)
% CORE_LOOP_HESS  Unified H2 mass/pressure/compressor/electrochem integrator.

time = time(:); P_el = P_el(:); P_fc = P_fc(:);
N   = numel(time);
dt  = time(2) - time(1);                     % hours (assumed constant, >0)

% ---------- preallocate ----------
n_H2        = zeros(N,1);    % mol
mass_H2     = zeros(N,1);    % kg
m_dot_el    = zeros(N,1);    % kg/h
m_dot_fc    = p.init.m_dot_fc(:);   % kg/h
P_tank      = zeros(N,1);    % Pa
P_tank_bar  = zeros(N,1);    % bar
P_comp      = zeros(N,1);    % kW

i_model     = zeros(N,1);    % A/cm^2
V           = zeros(N,1);    % V
eta_F       = zeros(N,1);    % -
eta_el      = zeros(N,1);    % LHV overall
P_el_remain = zeros(N,1);    % kW (power to EL after compressor)


% ---------- pressure-safe caps (ideal gas) ----------
m_cap_press = (p.P_outlet * p.tank.V / (p.R * p.tank.T)) * p.M_H2;  % kg
m_min = 0.05*m_cap_press; 
m_max = 0.95*m_cap_press;

% ---------- initial tank ----------
mass_H2(1) = min(max(p.init.mass_H2(1), m_min), m_max);
n_H2(1)    = mass_H2(1)/p.M_H2;
P_tank(1)  = n_H2(1)*p.R*p.tank.T/p.tank.V;
P_tank_bar(1) = P_tank(1)/1e5;


% ---------- main loop ----------
for k = 1:N-1
  % a) fuel cell draw (cap by available mass this step)
  if P_fc(k) > 0
    m_req    = P_fc(k)/(p.eta_fc*p.LHV_H2);     % kg/h
    m_cap_fc = max(0, (mass_H2(k)-m_min)/dt);   % kg/h
    m_dot_fc(k) = min(m_req, m_cap_fc);
  else
    m_dot_fc(k) = 0;
  end

  % b) available PV (kW) and pressure ratio (Pa/Pa)
  P_avail = max(0, P_el(k));                    % kW
  PR = (p.P_inlet>0) * (P_tank(k)/p.P_inlet); if PR<=1, PR=1; end

  % c) compressor + electrolyzer
 if p.enable.electrochem
    % single solve: P_el(i)+P_comp(i)=P_avail
    balance = @(i) ...
       (p.el.N*i*p.el.A*cell_voltage(i,p))/1e3 ...
     + ( (p.el.N*i*p.el.A*faraday_eff(i,p))/(2*p.F) * p.Cp_H2*p.tank.T ...
         * max(PR^((p.gamma-1)/p.gamma)-1,0) / p.eta_C )/1e3 ...
     - P_avail;

    if P_avail <= 0
      i_model(k) = 0;
    else
      % quick bound check at design current
      Pel_id   = (p.el.N*p.el.i*p.el.A*cell_voltage(p.el.i,p))/1e3;
      nH2_id   = (p.el.N*p.el.i*p.el.A*faraday_eff(p.el.i,p))/(2*p.F);
      Pcmp_id  = (nH2_id*p.Cp_H2*p.tank.T*max(PR^((p.gamma-1)/p.gamma)-1,0)/p.eta_C)/1e3;
      if Pel_id + Pcmp_id <= P_avail
        i_model(k) = p.el.i;
      else
        i_model(k) = fzero(@(i) balance(i), [0, p.el.i]);
      end
    end

    if i_model(k) <= 0
      % off
      V(k) = 0; eta_F(k)=0; m_dot_el(k)=0; P_comp(k)=0; P_el_remain(k)=0; eta_el(k)=0;
    else
      % on
      V(k)           = cell_voltage(i_model(k), p);
      eta_F(k)       = faraday_eff(i_model(k), p);
      nH2_s          = (p.el.N*i_model(k)*p.el.A*eta_F(k))/(2*p.F);     % mol/s
      P_el_remain(k) = (p.el.N*i_model(k)*p.el.A*V(k))/1e3;             % kW
      P_comp(k)      = ( nH2_s*p.Cp_H2*p.tank.T*max(PR^((p.gamma-1)/p.gamma)-1,0)/p.eta_C )/1e3;
      m_dot_el(k)    = nH2_s * p.M_H2 * 3600;                              % kg/h
      eta_el(k)      = (P_el_remain(k)>0) * (m_dot_el(k)*p.LHV_H2 / P_el_remain(k));
    end

  else
  % electrochem OFF (M1/M2)
 if P_avail <= 0
  m_dot_el(k)=0; P_comp(k)=0; V(k)=0; eta_F(k)=0; eta_el(k)=0; P_el_remain(k)=0;
 else
end
  if p.enable.compressor
    % shared power balance on m_dot_el (same compressor math as M3)
    % EL power:   m*(LHV/eta_el)
    % COMP power: ((m/M)/3600)*Cp*T*(PR^((gamma-1)/gamma)-1)/eta_C / 1e3
    m_dot_el(k) = fzero(@(m) ...
        ( m * (p.LHV_H2 / p.eta_el) ) + ...
        ( ((m/p.M_H2)/3600) * p.Cp_H2 * p.tank.T * max(PR^((p.gamma-1)/p.gamma)-1,0) / p.eta_C )/1e3 ...
        - P_avail, ...
        [0, P_avail/(p.LHV_H2/p.eta_el)]);

    nH2_s          = (m_dot_el(k)/p.M_H2)/3600;  % mol/s
    P_comp(k)      = ( nH2_s*p.Cp_H2*p.tank.T*max(PR^((p.gamma-1)/p.gamma)-1,0)/p.eta_C )/1e3;
    P_el_remain(k) = max(P_avail - P_comp(k), 0);   % kW to EL
    eta_el(k)      = (P_el_remain(k)>0) * p.eta_el;
    V(k)=0; eta_F(k)=0; i_model(k)=0;

  else
    % compressor OFF: all power to EL at fixed efficiency
    i_model(k)=0; V(k)=0; eta_F(k)=0;
    P_comp(k)=0;
    m_dot_el(k)    = p.eta_el * P_avail / p.LHV_H2;   % kg/h
    P_el_remain(k) = P_avail;                         % kW
    eta_el(k)      = (P_el_remain(k)>0) * p.eta_el;
  end
 end

  % --- enforce pressure window via FLOWS (not state) ---
  headroom       = max(0, m_max - mass_H2(k));       % kg available this step
  m_dot_el(k)    = min(m_dot_el(k), headroom / dt);  % kg/h
  if headroom == 0
      % full → no production/compression this step
      P_comp(k)=0; P_el_remain(k)=0; i_model(k)=0; V(k)=0; eta_F(k)=0; eta_el(k)=0;
  end

  % d) state update → k+1
  mass_H2(k+1) = mass_H2(k) + (m_dot_el(k) - m_dot_fc(k))*dt;
  mass_H2(k+1) = min(max(mass_H2(k+1), m_min), m_max);   % tiny safety
  n_H2(k+1)    = mass_H2(k+1)/p.M_H2;
  P_tank(k+1)  = n_H2(k+1)*p.R*p.tank.T/p.tank.V;
  P_tank_bar(k+1) = P_tank(k+1)/1e5;
end

% ---------- package ----------
out.time        = time;
out.P_el        = P_el;
out.P_fc        = P_fc;
out.mass_H2     = mass_H2;
out.m_dot_el    = m_dot_el;
out.m_dot_fc    = m_dot_fc;
out.P_comp_kW   = P_comp;           % keep external name
out.P_tank_bar  = P_tank_bar;
out.i_model     = i_model;
out.V           = V;
out.eta_F       = eta_F;
out.eta_el      = eta_el;
out.P_el_remain = P_el_remain;

end 