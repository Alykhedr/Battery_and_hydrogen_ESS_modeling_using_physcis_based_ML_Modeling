% pinn_soc_multigamma_ascii.m
% Trains a 2-input PINN [t, g_val]→SOC using your trapezoidal profile

clear; clc; close all;

%% 1) Cell & model parameters
Q_N      = 1.1;           % Ah
I_N      = 1.0;           % A
T_nom    = 25;            % °C
T_amb    = 22;            % °C
n        = -0.03;         % efficiency exponent
dt       = 1;             % hr
t_end    = 48;            % hr

g_vals = [1.00, 0.96, 0.94, 0.92];   % SOH fractions

%% 2) Build your trapezoidal current profile
T_vec   = (0:dt:t_end)';                % hours
N       = numel(T_vec);
pattern = [1, 0, -0.5, 0];              % [1h discharge,1h rest,1h charge,1h rest]
I_C     = repmat(pattern,1,ceil(N/numel(pattern)));
I_C     = I_C(1:N);
I_vec   = I_C * I_N;                    % current [A]

%% 3) PINN architecture: inputs [t_norm; g_val]
layers = [
    featureInputLayer(2,"Normalization","none","Name","in")
    fullyConnectedLayer(20,"Name","fc1")
    tanhLayer("Name","tanh1")
    fullyConnectedLayer(20,"Name","fc2")
    tanhLayer("Name","tanh2")
    fullyConnectedLayer(1,"Name","fc3")
    sigmoidLayer("Name","sig")   % SOC in (0,1)
];
dlnet = dlnetwork(layerGraph(layers));

%% 4) Training hyperparams
numEpochs = 3000;
learnRate = 1e-3;
Ncoll     = 1000;
lambdaIC  = 100;

%% 5) Training loop
avgG = []; avgG2 = [];
for ep = 1:numEpochs

    % a) Collocation batch (rows)
    t_coll = rand(1,Ncoll) * t_end;               % 1×Ncoll
    idx    = randi(numel(g_vals),1,Ncoll);
    g_coll = g_vals(idx);                         % 1×Ncoll

    % normalize time to [-1,1]
    tN = 2*(t_coll - t_end/2)/t_end;               % 1×Ncoll

    % assemble input
    Xcoll = dlarray([tN; g_coll],"CB");           % 2×Ncoll

    % b) IC batch at t=0→tN=-1, for all g_vals
    tIC = -1;                                     % 1×1
    GIC = g_vals;                                 % 1×numG
    Xic = dlarray([repmat(tIC,1,numel(GIC)); GIC],"CB");

    % c) compute loss + grads
    [loss, grads] = dlfeval(@modelLoss, dlnet, Xcoll, Xic, ...
                           Q_N, I_N, n, T_nom, T_amb, ...
                           T_vec, I_vec, t_end, lambdaIC);

    % d) Adam step
    [dlnet, avgG, avgG2] = adamupdate(dlnet, grads, avgG, avgG2, ep, learnRate);

    % e) display
    if mod(ep,500)==0
        fprintf("Epoch %4d | Loss = %.3e\n", ep, double(gather(extractdata(loss))));
    end
end

%% 6) Plot PINN vs True ODE for each g_val
figure; hold on; grid on;
tPlot = linspace(0,t_end,200);
tNplt = 2*(tPlot - t_end/2)/t_end;

for j = 1:numel(g_vals)
    g_val = g_vals(j);

    % PINN prediction
    Xp   = dlarray([tNplt; repmat(g_val,1,numel(tPlot))],"CB");
    SocP = forward(dlnet, Xp);
    SocP = gather(extractdata(SocP))';

    % True ODE (forward Euler)
    SocT = zeros(size(tPlot));
    SocT(1) = 1;
    for k = 2:length(tPlot)
        I_k  = interp1(T_vec, I_vec, tPlot(k));
        eta  = (abs(I_k)/I_N)^n;
        K_T  = 1 + 0.008*(T_nom - T_amb);
        dSoc = (K_T * eta * I_k * dt) / (g_val * Q_N);
        SocT(k) = max(0,min(1,SocT(k-1) - dSoc));
    end

    plot(tPlot, SocP,   '-', 'LineWidth',1.5, 'DisplayName',sprintf('PINN g=%.2f',g_val));
    plot(tPlot, SocT,  '--', 'LineWidth',1.5, 'DisplayName',sprintf('True g=%.2f',g_val));
end

xlabel('Time [h]');
ylabel('SOC (0–1)');
title('PINN vs True SOC for Multiple SOH (g)');
legend('Location','eastoutside');


%% --- modelLoss subfunction ---
function [loss, grads] = modelLoss(dlnet, Xc, Xic, ...
                                   Q_N, I_N, n, Tnom, Tamb, ...
                                   Tvec, Ivec, t_end, λIC)

    % 1) Predict SOC at collocation
    SOCc = forward(dlnet, Xc);                % 1×Ncoll

    % 2) dSOC/d(t_norm)
    dSdtn = dlgradient(sum(SOCc,"all"), Xc);  % 2×Ncoll
    dSdtn = dSdtn(1,:);                       % wrt time

    % chain rule → dSOC/dt
    dSdt = dSdtn * (2/t_end);

    % 3) recover real t & g_val
    t_norm = extractdata(Xc(1,:))';           % Ncoll×1
    t_real = 0.5*(t_norm + 1)*t_end;
    g_coll = extractdata(Xc(2,:))';           % Ncoll×1

    % 4) I and eta
    Icol = interp1(Tvec, Ivec, t_real);
    eta  = (abs(Icol)/I_N).^n;

    % 5) thermal factor
    K_T  = 1 + 0.008*(Tnom - Tamb);

    % 6) ODE residual
    R = dSdt + (K_T .* eta .* Icol)' ./ (g_coll' * Q_N);

    % 7) IC loss: SOC(0)=1 for all g_vals
    Sic = forward(dlnet, Xic);                % 1×numG
    lossIC = mean((Sic - 1).^2,"all");

    % 8) residual loss
    lossRes = mean(R.^2,"all");

    % total
    loss  = lossRes + λIC * lossIC;
    grads = dlgradient(loss, dlnet.Learnables);
end
