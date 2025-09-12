%% File: plot_bess.m
function plot_bess(time, P_pv, P_load, SOC, SoH) %#ok<INUSD>
% Two stacked plots:
% 1) PV vs Demand (time)
% 2) SOC (time)

    figure('Color','w','Position',[100 100 1100 650]);

    % ---- (1) PV & Demand (time) ----
    ax1 = subplot(2,1,1);
    plot(ax1, time, P_pv,'b'); hold(ax1,'on');
    plot(ax1, time, P_load,'r');
    grid(ax1,'on'); ylabel(ax1,'Power (kW)');
    title(ax1,'PV Generation & Demand');
    legend(ax1, 'PV gen','Demand','Location','northwest');

    % ---- (2) SOC (time) ----
    ax2 = subplot(2,1,2);
    plot(ax2, time, 100*SOC, 'LineWidth',0.9);
    grid(ax2,'on'); ylabel(ax2,'SOC (%)'); ylim(ax2,[0 100]);
    title(ax2,'Battery State of Charge');

    % Keep time axes aligned
    linkaxes([ax1 ax2],'x');

    % ----- Tick formatting depends on time type -----
    if isdatetime(time)
        xtickformat(ax1,'MMM yyyy');
        xtickformat(ax2,'MMM yyyy');
        xlabel(ax2,'Date');
    elseif isduration(time)
        % duration axis (e.g., hours/minutes)
        xtickformat(ax1,'hh:mm');
        xtickformat(ax2,'hh:mm');
        xlabel(ax2,'Time');
    else
        % numeric hours
        xtickformat(ax1,'%.0f');
        xtickformat(ax2,'%.0f');
        xlabel(ax2,'Time (h)');
    end
end
