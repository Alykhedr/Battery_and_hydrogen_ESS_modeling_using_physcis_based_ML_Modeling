function Ua = Ua_graphite(soc, map)
% UA_GRAPHITE  Graphite anode OCV vs SOC using Safari et al. form with linear SOC→xa map.
% map.xa0 and map.xa100 can be provided (calibrated). Defaults provided if missing.

if nargin < 2 || ~isstruct(map)
    map.xa0   = 0.01;
    map.xa100 = 0.78992288; % calibrated so Ua(50%) ~ 0.123 V
end

% Clamp SOC
soc = max(0, min(1, soc));

% Map SOC → anode lithiation xa
xa = map.xa0 + soc*(map.xa100 - map.xa0);

% Safari et al. anode OCV (V)
Ua = 0.6379 ...
   + 0.5416*exp(-305.5309*xa) ...
   + 0.0440*tanh(-(xa - 0.1958)/0.1088) ...
   - 0.1978*tanh((xa - 1.0571)/0.0854) ...
   - 0.6875*tanh((xa + 0.0117)/0.0529) ...
   - 0.0175*tanh((xa - 0.5692)/0.0875);
end
