function out = core_loop(time, Ppv, Pload, p)
% CORE_LOOP  Unified SOC/SOH/EFC integrator for B1/B2/B3 via param gating.

time = time(:); Ppv = Ppv(:); Pload = Pload(:);
N = numel(time);

% Numeric hours from start
if isdatetime(time)
    time_h = hours(time - time(1));
elseif isduration(time)
    time_h = hours(time);
else
    time_h = time - time(1); % assume hours
end

% Prealloc
Q      = zeros(N,1);     % Ah stored charge
Ca     = zeros(N,1);     % Ah available capacity
SOC    = zeros(N,1);     % 0..1
SOH    = ones(N,1);      % 0..1
P_bess = zeros(N-1,1);   % kW (+dis / -chg)
EFC    = zeros(N,1);     % discharge-only full equivalent cycles

% Loss trackers (percent points)
cyc_pct = zeros(N,1);
cal_pct = zeros(N,1);

% ---- Accumulators must exist regardless of gating ----
Qloss_cyc_pct = 0;   % cumulative %-points from cycling
Qloss_cal_pct = 0;   % cumulative %-points from calendar
EFC_prev      = 0;   % for Δ(EFC^a)

% Initial conditions
Ca(1)  = p.Q_nom;              % Ah
Q(1)   = p.SOC0 * Ca(1);
SOC(1) = Q(1)/Ca(1);

% --- B1 energy state (only if BOTH cycle & calendar are OFF)
isB1 = ~p.modes.do_cycle && ~p.modes.do_calendar;
if isB1
    E    = zeros(N,1);                  % kWh
    Enom = p.V_nom * p.Q_nom / 1000;    % kWh
    E(1) = p.SOC0 * Enom;
    Emin = p.SOC_min * Enom;
    Emax = p.SOC_max * Enom;
end


% Precompute Montes factors if cycle enabled
if p.modes.do_cycle
    T_K     = p.T_degC + 273.15;
    fT_cyc  = exp( p.kT * (T_K - p.Tref_K) / T_K );
    DOD_pct = 100 * (p.SOC_max - p.SOC_min);
    fD_cyc  = exp( p.kDODc * DOD_pct );
    % Discharge-weighted mSOC proxy
    soc_wsum = 0.0; ah_wsum = 0.0; mSOC_eff = SOC(1);
end

% Calendar accumulators if enabled
if p.modes.do_calendar
    tcal_h = 0;
    % cache constants
    fTcal_const = exp( -p.deg.Ea_cal_Jmol/p.R_gas * (1/(p.T_degC+273.15) - 1/p.deg.Tref_K) );
end

% ---- Main loop
for k = 1:N-1
    h = max(time_h(k+1) - time_h(k), 0);
    Pnet = Pload(k) - Ppv(k);                 % +dis / -chg


        % ===== B1: energy-balance path (no cycle, no calendar) =====
    if isB1
        % Ask, capped by hardware
        Pask = min(max(Pnet, -p.P_ch_max), p.P_dis_max);  % kW
    
        % SOC-window headroom (kW) for this step
        if Pask >= 0
            % discharge
            Plim_soc = max(0, (E(k) - Emin) * p.eta_dis / h);
            P        = min(Pask, Plim_soc);
        else
            % charge
            Plim_soc = max(0, (Emax - E(k)) / (p.eta_ch * h));
            P        = max(Pask, -Plim_soc);
        end
    
        % Commit power & update energy (only place with η in state update)
        P_bess(k) = P;
        E(k+1) = min(max( E(k) + ( max(-P,0)*p.eta_ch - max(P,0)/p.eta_dis )*h , Emin), Emax);
    
        % Derive SOC and keep Q/Ca coherent for plots
        SOC(k+1) = E(k+1)/Enom;
        Ca(k+1)  = p.Q_nom;
        Q(k+1)   = SOC(k+1) * Ca(k+1);
    
        % B1 has no cycle/calendar or SOH updates
        continue
    end
    % ===== end B1 branch =====

    % Power limit (kW) and requested current (A)
    Plim = min(max(Pnet, -p.P_ch_max), p.P_dis_max);
    Ireq = 1000*Plim / p.V_nom;               % A (+dis, -chg)

    % Requested split and ΔQ (Ah)
    Ich_req  = max(-Ireq,0);
    Idis_req = max( Ireq,0);
    dQ_req   = h*( Ich_req - Idis_req );

    % SOC-window clamp on Q
    Qmin = p.SOC_min*Ca(k);  Qmax = p.SOC_max*Ca(k);
    Qn   = min(max(Q(k)+dQ_req, Qmin), Qmax);
    dQa  = Qn - Q(k);

    % Actual currents from clamped ΔQ
    if h > 0
        if dQa >= 0
            Ich_a =  dQa/h;  Idis_a = 0;
        else
            Ich_a =  0;      Idis_a = (-dQa)/h;
        end
    else
        Ich_a = 0; Idis_a = 0;
    end

    % Net power (efficiency only)
    P_bess(k) = (p.V_nom*Idis_a)/1000 * p.eta_dis ...
              - (p.V_nom*Ich_a )/1000 / p.eta_ch;

    % Update stored charge
    Q(k+1) = Qn;

    % ------------- CYCLE AGING (Montes) -------------
    if p.modes.do_cycle
        % Discharge-only EFC
        EFC(k+1) = EFC(k) + (h * Idis_a) / p.Q_nom;

        % Per-step C-rates
        Cdis = Idis_a / p.Q_nom;
        Cch  = Ich_a  / p.Q_nom;

        % Discharge-weighted mSOC
        if Idis_a > 0
            dAh_dis  = h * Idis_a;
            soc_wsum = soc_wsum + SOC(k) * dAh_dis;
            ah_wsum  = ah_wsum  + dAh_dis;
            mSOC_eff = soc_wsum / max(ah_wsum, eps);
        end
        mSOC = mSOC_eff;

        % δ_cyc (pp per ΔEFC^a)
        delta = cycle_aging(Cch, Cdis, mSOC, p, fT_cyc, fD_cyc);

        % Increment consistent with Cf = δ * EFC^a
        EFC_now        = EFC(k+1);
        dQloss_cyc_pp  = delta * ( EFC_now^p.a_montes - EFC_prev^p.a_montes );
        Qloss_cyc_pct  = Qloss_cyc_pct + dQloss_cyc_pp;
        EFC_prev       = EFC_now;
    else
        % keep EFC flat and cycle loss at zero when disabled
        EFC(k+1) = EFC(k);
    end

    % ------------- CALENDAR AGING (Naumann/Schimpe) -------------
    if p.modes.do_calendar
        [dQcal_pp, ~, ~] = calendar_aging(SOC(k), h, tcal_h, p, fTcal_const);
        Qloss_cal_pct = Qloss_cal_pct + dQcal_pp;
        tcal_h        = tcal_h + h;
    end

    % Trackers (safe: both accumulators always exist)
    cyc_pct(k+1) = Qloss_cyc_pct;
    cal_pct(k+1) = Qloss_cal_pct;

    % SOH/Ca update & SOC continuity
    loss_total = cyc_pct(k+1) + cal_pct(k+1);
    SOH(k+1) = max(p.deg.SOH_min, 1 - loss_total/100);
    Ca_old   = Ca(k);
    Ca(k+1)  = SOH(k+1) * p.Q_nom;

    SOC_cont = Q(k+1)/max(Ca_old, eps);
    Qtarget  = SOC_cont * Ca(k+1);
    Q(k+1)   = min(max(Qtarget, p.SOC_min*Ca(k+1)), p.SOC_max*Ca(k+1));
    SOC(k+1) = Q(k+1)/Ca(k+1);
end

% ---- Package outputs
out.time    = time;
out.SOC     = SOC;
out.SOH     = SOH;
out.EFC     = EFC;
out.P_bess  = P_bess;
out.loss.cycle    = cyc_pct;
out.loss.calendar = cal_pct;
out.meta.model = p.meta.model;
out.meta.params = p;
end
