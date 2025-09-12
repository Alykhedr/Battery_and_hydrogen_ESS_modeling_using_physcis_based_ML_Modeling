%% run_models.m
% Wrapper to run Level1–3 HESS models, compute comparison metrics, display table, and plot

  clear; clc;

  % Define metric and model names
  varNames   = {'InitMass','FinalMass','TotProd','TotCons','PeakProd', ...
                'Comp_E_kWh','PV2H2_pct','RTE_pct'};
  modelNames = {'Level1','Level2','Level3'};
  nModels    = numel(modelNames);
  metrics    = zeros(nModels, numel(varNames));
  tankMass   = cell(1,nModels);
  timeVec    = [];

  LHV_H2 = 33.33;  % kWh/kg

  % Loop over models 1–3
  for m = 1:nModels
    fprintf('Running %s...\n', modelNames{m});
    switch m
      case 1
        [time,P_el,P_fc,m_dot_el,m_dot_fc,P_comp_kW,~,mass_H2] = run_HESS1();
      case 2
        [time,P_el,P_fc,m_dot_el,m_dot_fc,P_comp_kW,~,mass_H2] = run_HESS2();
      case 3
        [time,P_el,P_fc,m_dot_el,m_dot_fc,P_comp_kW,eta_F,mass_H2] = run_HESS3();
    end

    % store time & tank mass
    if isempty(timeVec)
      timeVec = time;
    end
    tankMass{m} = mass_H2;

    % Compute metrics
    InitMass   = mass_H2(1);
    FinalMass  = mass_H2(end);
    TotProd    = trapz(time, m_dot_el);
    TotCons    = trapz(time, m_dot_fc);
    PeakProd   = max(m_dot_el);
    PV_E_kWh   = trapz(time, P_el);
    Comp_E_kWh = trapz(time, P_comp_kW);
    PV2H2_pct  = TotProd * LHV_H2 / PV_E_kWh * 100;
    FC_Energy  = trapz(time, P_fc);
    RTE_pct    = FC_Energy / (PV_E_kWh + Comp_E_kWh) * 100;

    metrics(m,:) = [InitMass,FinalMass,TotProd,TotCons,PeakProd,...
                    Comp_E_kWh,PV2H2_pct,RTE_pct];
  end

  % Display table
  T = array2table(metrics,'RowNames',modelNames,'VariableNames',varNames);
  disp(T);

  % Plot 1: Tank mass comparison (Models 1–3)
  figure('Name','Tank Mass Comparison','Units','normalized','Position',[.1 .1 .6 .4]);
  hold on;
    for m=1:nModels
      plot(timeVec, tankMass{m}, 'LineWidth',1.5);
    end
  hold off;
  xlabel('Time (hr)'); ylabel('H_2 Tank Mass (kg)');
  title('H_2 Tank Mass: Models 1–3');
  legend(modelNames,'Location','best'); grid on;

    % Plot 2: Bar chart for PV→H₂ conversion and RTE (all three models)
  figure('Name','Efficiency Comparison','Units','normalized','Position',[.1 .1 .6 .4]);
  h = bar(metrics(:,[7,8]), 'grouped');   % column 8 = PV2H2_pct, column 9 = RTE_pct
  set(gca,'XTick',1:nModels,'XTickLabel',modelNames);
  xlabel('HESS Model');
  ylabel('Efficiency (%)');
  title('PV→H₂ Conversion vs Round-Trip Efficiency');
  legend(h, {'PV→H₂ (%)','RTE (%)'}, ...
         'Location','northoutside','Orientation','horizontal');
  grid on;



%% Local functions for each model
function [time,P_el,P_fc,m_dot_el,m_dot_fc,P_comp_kW,eta_F,mass_H2] = run_HESS1()
  [time,P_el,P_fc] = load_profiles();
  p = params(time,1); p.enable.compressor = false;
  N = numel(time); dt = time(2)-time(1);
  m_dot_el   = p.init.m_dot_el;
  m_dot_fc   = p.init.m_dot_fc;
  mass_H2    = p.init.mass_H2;
  eta_F      = ones(1,N);
  P_comp_kW  = zeros(1,N);
  for k = 2:N
    if P_el(k) > 0
      m_dot_el(k) = p.eta_el * P_el(k) / p.LHV_H2;
    end
    if P_fc(k) > 0
      m_dot_fc(k) = P_fc(k) / (p.eta_fc * p.LHV_H2);
    end
    mass_H2(k) = mass_H2(k-1) + (m_dot_el(k) - m_dot_fc(k)) * dt;
  end
end

function [time,P_el,P_fc,m_dot_el,m_dot_fc,P_comp_kW,eta_F,mass_H2] = run_HESS2()
  [time,P_el,P_fc] = load_profiles(); p = params(time,2);
  N = numel(time); dt = time(2)-time(1);
  n_H2      = zeros(1,N);
  P_comp    = zeros(1,N);
  m_dot_el  = p.init.m_dot_el;
  m_dot_fc  = p.init.m_dot_fc;
  eta_F     = ones(1,N);
  n_H2(1)   = p.init.mass_H2(1) / p.M_H2;
  for k = 2:N
    if P_fc(k)>0
      m_dot_fc(k) = P_fc(k)/(p.eta_fc*p.LHV_H2);
    else
      m_dot_fc(k) = 0;
    end
    n_H2(k) = n_H2(k-1)+(m_dot_el(k-1)-m_dot_fc(k))*dt/p.M_H2;
    PR = p.P_outlet/max((n_H2(k)*p.R*p.tank.T/p.tank.V),p.P_inlet);
    mdot_mol = (p.eta_el*P_el(k)/p.LHV_H2)/p.M_H2/3600;
    W = mdot_mol*p.Cp_H2*p.tank.T*(PR^((p.gamma-1)/p.gamma)-1)/p.eta_C;
    P_comp(k)=max(0,W);
    used = min(P_comp(k)/1e3,P_el(k));
    m_dot_el(k)=p.eta_el*(P_el(k)-used)/p.LHV_H2;
  end
  mass_H2   = n_H2*p.M_H2;
  P_comp_kW = P_comp/1e3;
end

function [time,P_el,P_fc,m_dot_el,m_dot_fc,P_comp_kW,eta_F,mass_H2] = run_HESS3()
  [time,P_el,P_fc] = load_profiles(); p = params(time,3);
  N = numel(time); dt = time(2)-time(1);
  n_H2     = zeros(1,N);
  P_comp   = zeros(1,N);
  m_dot_el = p.init.m_dot_el;
  m_dot_fc = p.init.m_dot_fc;
  eta_F    = zeros(1,N);
  area_cm2 = p.el.A*1e4;
  n_H2(1)  = p.init.mass_H2(1)/p.M_H2;
  for k = 2:N
    if P_fc(k)>0
      m_dot_fc(k)=P_fc(k)/(p.eta_fc*p.LHV_H2);
    else
      m_dot_fc(k)=0;
    end
    P_av = max(0,P_el(k));
    Vd   = cell_voltage(p.el.i,p);
    P_req= p.el.N*p.el.i*area_cm2*Vd/1e3;
    frac = min(1,P_av/P_req);
    im   = frac*p.el.i;
    I_tot= p.el.N*im*area_cm2;
    ef   = faraday_eff(im,p);
    PR   = p.P_outlet/max((n_H2(k-1)*p.R*p.tank.T/p.tank.V),p.P_inlet);
    Wc   = I_tot*ef/(2*p.F)*p.Cp_H2*p.tank.T*(PR^((p.gamma-1)/p.gamma)-1)/p.eta_C;
    P_comp(k) = min(P_av,Wc/1e3)*1e3;
    P_rem = P_av - P_comp(k)/1e3;
    if P_rem>0
      fun = @(i)p.el.N*i*area_cm2.*cell_voltage(i,p)/1e3 - P_rem;
      ia  = fzero(fun,[0,p.el.i]);
      ia  = max(0,min(p.el.i,ia));
      I_tot= p.el.N*ia*area_cm2;
      ef   = faraday_eff(ia,p);
      m_dot_el(k) = I_tot*ef/(2*p.F)*p.M_H2*3600;
      eta_F(k)    = ef;
    end
    n_H2(k)= min(max(n_H2(k-1)+(m_dot_el(k)-m_dot_fc(k))*dt/p.M_H2,0),p.max_mass_H2/p.M_H2);
  end
  mass_H2   = n_H2*p.M_H2;
  P_comp_kW = P_comp/1e3;
end
