% verify_model1.m — Minimal, readable checks for B1 (energy balance only)
% Tests:
%  A) Zero net power → SOC constant
%  B) Constant discharge (below limits)
%  C) Constant charge (below limits)  [PV>0, Load=0]
%  D) Bounds clamp at SOC walls (SOC_max / SOC_min)
%  E) Caps: hardware limit vs SOC headroom
%  F) Energy residual ~ 0 (conservation)

clear; clc;

% ---------- Params ----------
p = params('B1');                       % cycle/calendar OFF
Enom_kWh = (p.V_nom * p.Q_nom)/1000;    % nominal energy
tol.soc_abs   = 1e-6;
tol.pow_abs   = 1e-9;
tol.conv_soc  = 5e-4;

say = @(varargin) fprintf('%s\n', sprintf(varargin{:}));

results = struct('name',{},'pass',{},'info',{});

%% A) Zero net power
name = 'A Zero net power';
SOC0 = 0.50; T=10; dt=1; t=(0:dt:T)';
Ppv   = 100*ones(size(t));
Pload = 100*ones(size(t));              % Pnet = 0
pA = p; pA.SOC0 = SOC0;
out = core_loop(t, Ppv, Pload, pA);
pass = max(abs(out.SOC - SOC0)) <= tol.soc_abs;
results(end+1) = pack(name, pass, struct('maxAbsSOCerr', max(abs(out.SOC - SOC0))));

%% B) Constant discharge (below limits)
name = 'B Const discharge';
SOC0 = 0.8; Pdis=100; T=4; dt=1; t=(0:dt:T)';
Ppv   = zeros(size(t));
Pload = Pdis*ones(size(t));             % Pnet = +Pdis → discharge
pB = p; pB.SOC0 = SOC0;
out = core_loop(t, Ppv, Pload, pB);
SOC_exp = SOC0 - (Pdis/p.eta_dis)*T/Enom_kWh;
pass = abs(out.SOC(end)-SOC_exp) <= tol.soc_abs;
results(end+1) = pack(name, pass, struct('SOC_end',out.SOC(end),'SOC_expected',SOC_exp));

%% C) Constant charge (below limits) [PV>0, Load=0]
name = 'C Const charge';
SOC0 = 0.2; Pchg=80; T=4; dt=1; t=(0:dt:T)';
Ppv   = Pchg*ones(size(t));             % PV supplies charge
Pload = zeros(size(t));                 % Pnet = -Pchg → charge
pC = p; pC.SOC0 = SOC0;
out = core_loop(t, Ppv, Pload, pC);
SOC_exp = SOC0 + (Pchg*p.eta_ch)*T/Enom_kWh;
pass = abs(out.SOC(end)-SOC_exp) <= tol.soc_abs;
results(end+1) = pack(name, pass, struct('SOC_end',out.SOC(end),'SOC_expected',SOC_exp));

%% D) Bounds clamp at SOC walls
name = 'D Bounds clamp';
dt=1; T=3; t=(0:dt:T)';
% E1: at SOC_max try to charge → no motion, P≈0
pE1 = p; pE1.SOC0 = p.SOC_max;
out1 = core_loop(t, Pchg_profile(200,t), zeros(size(t)), pE1);
ok1 = max(abs(out1.P_bess))<=tol.pow_abs && (max(out1.SOC)-min(out1.SOC))<=tol.soc_abs;
% E2: at SOC_min try to discharge → no motion, P≈0
pE2 = p; pE2.SOC0 = p.SOC_min;
out2 = core_loop(t, zeros(size(t)), Pdis_profile(200,t), pE2);
ok2 = max(abs(out2.P_bess))<=tol.pow_abs && (max(out2.SOC)-min(out2.SOC))<=tol.soc_abs;
pass = ok1 && ok2;
results(end+1) = pack(name, pass, struct());

%% E) Caps: hardware vs headroom
% E1: hardware bind (use short dt so SOC headroom is wide)
name = 'E Caps';
dt=0.1; T=0.1; t=(0:dt:T)'; SOC0=0.5;
% Charge: ask huge PV; expect P ≈ -p.P_ch_max on 1st step
pF = p; pF.SOC0 = SOC0;
Ppv = 10*p.P_ch_max*ones(size(t)); Pload = zeros(size(t));
outF = core_loop(t, Ppv, Pload, pF);
hw_ok_ch = abs(outF.P_bess(1) + p.P_ch_max) <= 1e-6;
% E2: headroom bind (use long dt near SOC_max)
dt=1; T=1; t=(0:dt:T)'; SOC0=p.SOC_max-1e-3;
pFH = p; pFH.SOC0 = SOC0;
Ppv = 10*p.P_ch_max*ones(size(t)); Pload = zeros(size(t));
outH = core_loop(t, Ppv, Pload, pFH);
E    = SOC0*Enom_kWh;
Emax = p.SOC_max*Enom_kWh;
P_soc = -(Emax - E)/(p.eta_ch*dt);  % expected headroom cap (negative)
head_ok = abs(outH.P_bess(1) - P_soc) <= 1e-6;
results(end+1) = pack(name, hw_ok_ch && head_ok, struct('P_bess_hw',outF.P_bess(1),'P_bess_hr',outH.P_bess(1),'P_soc_expect',P_soc));



%% F) Energy residual ~ 0 (conservation)
SOC0 = 0.73; Pdis=120; T=5; dt=0.25;
t = (0:dt:T)'; Ppv = zeros(size(t)); Pload = Pdis*ones(size(t));
pK = p; pK.SOC0 = SOC0;
outK = core_loop(t, Ppv, Pload, pK);
Enom = (p.V_nom*p.Q_nom)/1000;
dSOC = outK.SOC(end) - SOC0;

% Reconstruct energy flow using the model’s own sign/efficiency convention
E_flow = 0;
for k=1:numel(t)-1
    dth = t(k+1)-t(k);
    Pk  = outK.P_bess(k);  % +dis, -chg
    E_flow = E_flow + ( max(-Pk,0)*p.eta_ch - max(Pk,0)/p.eta_dis ) * dth;
end
residual = Enom*dSOC - E_flow;
passK = abs(residual) <= 1e-9;
results(end+1) = pack('F energy residual', passK, struct('residual',residual));

% ---------- Summary ----------
say('=== B1 Verification Summary (simple) ===');
for i=1:numel(results)
    tag = tern(results(i).pass,'PASS','FAIL');
    say('[%s] %s', tag, results(i).name);
end

% ---------- Helpers ----------
function r = pack(name,pass,info)
r.name=name; r.pass=logical(pass); if nargin<3, info=struct(); end; r.info=info;
end
function s = tern(c,a,b), if c, s=a; else, s=b; end
end
function P = Pchg_profile(Pchg,t), P = Pchg*ones(size(t)); end
function P = Pdis_profile(Pdis,t), P = Pdis*ones(size(t)); end

