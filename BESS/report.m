function report(out)
% REPORT  Print concise summary with losses and EFC.

time = out.time;
if isdatetime(time)
    total_h = hours(time(end) - time(1));
elseif isduration(time)
    total_h = hours(time(end));
else
    total_h = time(end) - time(1); % hours
end
span_days = total_h/24;

cyc = out.loss.cycle(end);
cal = out.loss.calendar(end);
tot = cyc + cal;

fprintf('Final SOH = %.2f%%  |  EFC = %.0f\n', 100*out.SOH(end), out.EFC(end));
fprintf('--- SUMMARY ---\n');
fprintf('Span: %.1f days (%.2f years)\n', span_days, span_days/365);
fprintf('Cycle loss    : %.2f%%\n', cyc);
fprintf('Calendar loss : %.2f%%\n', cal);
fprintf('Total loss    : %.2f%%\n', tot);
fprintf('EFC total     : %.0f (%.2f per day)\n', out.EFC(end), out.EFC(end)/max(span_days,eps));
end
