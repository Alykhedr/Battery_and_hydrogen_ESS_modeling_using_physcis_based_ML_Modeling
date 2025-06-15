function plot_hess(time, P_el, P_fc, mass_H2)
  figure;
  subplot(2,1,1)
    plot(time, P_el, '-b', time, P_fc, '-r');
    title('Power Profile');
    ylabel('Power [kW]'); legend('P_{el}','P_{fc}'); grid on;

  subplot(2,1,2)
    plot(time, mass_H2, '-k');
    title('H_2 Tank Level');
    ylabel('Mass [kg]'); xlabel('Time [hr]'); grid on;
end