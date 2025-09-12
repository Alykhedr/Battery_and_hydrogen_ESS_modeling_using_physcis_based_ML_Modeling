function [time, P_el, P_fc] = Loadprofileweek()
%LOAD_PROFILES  0–170 hr PV & FC demand matching the paper’s windows

  %% 1) Time axis
  dt   = 1;          
  T    = 167;        
  time = (0:dt:T)';  % column vector

  %% 2) PV injection P_el: single‐hump sine, 0 until t=5, peak=1000 kW
  P_el = 1000 * sin(2*pi*(time-5)/24);
  P_el(time < 5)    = 0;
  P_el(P_el  < 0)   = 0;

  %% 3) Fuel-cell demand P_fc: exact piecewise from your spec
  P_fc = zeros(size(time));

  P_fc(time>=  8  & time< 16) = 200;   %  8–16 h
  P_fc(time>= 16  & time< 32) =  50;   % 16–32 h
  P_fc(time>= 32  & time< 40) = 200;   % 32–40 h
  P_fc(time>= 40  & time< 41) = 125;   % 40–41 h
  P_fc(time>= 41  & time< 46) =   0;   % 41–46 h

  P_fc(time>= 46  & time< 56) =  50;   % 46–56 h
  P_fc(time>= 56  & time< 64) = 200;   % 56–64 h
  P_fc(time>= 64  & time< 65) = 150;   % 64–65 h
  P_fc(time>= 65  & time< 70) =   0;   % 65–70 h

  P_fc(time>= 70  & time< 80) =  50;   % 70–80 h
  P_fc(time>= 80  & time< 86) = 200;   % 80–86 h
  P_fc(time>= 86  & time< 88) = 125;   % 86–88 h
  P_fc(time>= 88  & time<104) =   0;   % 88–104 h

  P_fc(time>=104  & time<106) = 150;   %104–106 h
  P_fc(time>=106  & time<112) = 200;   %106–112 h
  P_fc(time>=112  & time<128) =   0;   %112–128 h

  P_fc(time>=128  & time<133) = 300;   %128–133 h
  P_fc(time>=133  & time<138) =   0;   %133–138 h
  P_fc(time>=138  & time<153) =  75;   %138–153 h
  P_fc(time>=153  & time<160) = 200;   %153–160 h

  P_fc(time>=160  & time<=170)=   0;   %160–170 h

end
