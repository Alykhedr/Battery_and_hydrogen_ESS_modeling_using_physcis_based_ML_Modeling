function plot_hess(time, P_el, P_fc, mass_H2, m_dot_el, m_dot_fc, varargin)
% PLOT_HESS  Plot HESS results with optional compressor & electrochem panels.
%   plot_hess(time,P_el,P_fc,mass_H2,m_dot_el,m_dot_fc)
%       -> 3-panel plot
%
%   plot_hess(...,P_comp_kW,P_tank_bar)
%       -> +1 panel for compressor & pressure
%
%   plot_hess(...,P_comp_kW,P_tank_bar,i_sol,V_cell,eta_F)
%       -> +2 panels for electrolyzer current & voltage/efficiency

  % parse optional inputs
  doComp = false;
  doElec = false;
  nOpt = numel(varargin);
  if nOpt == 2
    [P_comp_kW, P_tank_bar] = deal(varargin{:});
    doComp = true;
  elseif nOpt == 5
    [P_comp_kW, P_tank_bar, i_sol, V_cell, eta_F] = deal(varargin{:});
    doComp = true;
    doElec = true;
  end

  % total rows
  nrows = 3 + double(doComp) + double(doElec);

  figure('Units','normalized','Position',[.1 .1 .8 .8]);

  % Panel 1: Power Profile
  subplot(nrows,1,1);
  plot(time, P_el, '-b', time, P_fc, '-r','LineWidth',1.9);
  title('Power Profile'); ylabel('Power [kW]');
  legend('P_{el}','P_{fc}','Location','best'); grid on;
  xticks(0:10:time(end));

  % Panel 2: Tank Mass
  subplot(nrows,1,2);
  plot(time, mass_H2, '-k','LineWidth',1.9);
  title('H2 Tank Level'); ylabel('Mass [kg]');
  grid on; ylim([0 1.1*max(mass_H2)]);
  xticks(0:10:time(end));

  % Panel 3: H2 Flow Rates
  subplot(nrows,1,3);
  plot(time, m_dot_el, '-g', time, m_dot_fc, '-m','LineWidth',1.9);
  title('Hydrogen Flow Rates'); ylabel('kg/h'); xlabel('Time [hr]');
  legend('Production','Consumption','Location','best'); grid on;
  xticks(0:10:time(end));

  row = 4;
  % Panel 4: Compressor & Pressure
  if doComp
    subplot(nrows,1,row); row = row + 1;
    yyaxis left
      plot(time, P_comp_kW, '-m','LineWidth',1.9);
      ylabel('Compressor Power [kW]');
    yyaxis right
      plot(time, P_tank_bar, '--g','LineWidth',1.9);
      ylabel('Tank Pressure [bar]');
    title('Compressor & Tank Pressure');
    xlabel('Time [hr]'); legend('P_{comp}','P_{tank}','Location','best');
    grid on; xticks(0:10:time(end));
  end

  % Panels 5–6: Electrochemical variables
  if doElec
    % Panel 5: Current density
    subplot(nrows,1,row); row = row + 1;
    plot(time, i_sol, '-b','LineWidth',1.9);
    title('Electrolyzer Current Density');
    ylabel('Current [A/cm^2]'); xlabel('Time [hr]'); grid on;
    xticks(0:10:time(end));

    % Panel 6: Cell voltage & Faraday efficiency
    subplot(nrows,1,row);
    yyaxis left
      plot(time, V_cell, '-k','LineWidth',1.9);
      ylabel('Cell Voltage [V]');
    yyaxis right
      plot(time, eta_F, '--r','LineWidth',1.9);
      ylabel('Faraday Eff. [-]');
    title('Electrolyzer Voltage & Efficiency');
    xlabel('Time [hr]'); legend('Voltage','eta_F','Location','best');
    grid on; xticks(0:10:time(end));
  end

  % Always show full results table
  T = table( ...
    time(:), ...
    m_dot_el(:), ...
    m_dot_fc(:), ...
    mass_H2(:), ...
    'VariableNames', {'Hour','H2_prod','H2_cons','Tank_mass'} ...
  );
  disp(T);
end
