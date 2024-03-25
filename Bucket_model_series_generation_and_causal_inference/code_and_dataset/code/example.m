%% This script gives an example of  bucket model series generation and causal inference using the Liang-Kleeman information flow
%% It is written by Yuhui Zhao (yuhuizhao73@foxmail.com).

%% Parameters setting
num_steps = 10000; % number of time steps in one simulation
nsim = 100; % number of simulations
dt = 1; % time interval
IF_length = 1000; % sample length adopted for causal inference: 100, 300, 500, 1000, 1500, 2000, 2500, 3000
pn_SNR = 10000; % signal-noise-ratio of the process noise added: 2, 3, 4, 5, 10, 20, 50, 100, 10000
on_SNR = 10; % signal-noise-ratio of the observation noise added: 2, 3, 4, 5, 10, 20, 50, 100, 10000
on_switch = 1; % 1 (0) for construction with (without) observation noise
 
%% Deterministic simulation for the calculation of signal powers, with no 
%noise added 
[R, S, I, Q] = bucket_model_construction(num_steps, pn_SNR, 0, 0, 0);
power_R = var(R);
power_S = var(S);
power_I = var(I);
power_Q = var(Q);
 
%% Causal inference and inference evaluation
FasR = 0; % false rate
OmiR = 0; % omission rate
for isim = 1:nsim
    [R, S, I, Q] = bucket_model_construction(num_steps, pn_SNR, power_S, power_I, power_Q);
    % add observation noise
    ome_R = sqrt(power_R/on_SNR)*randn(num_steps, 1);
    ome_S = sqrt(power_S/on_SNR)*randn(num_steps, 1);
    ome_I = sqrt(power_I/on_SNR)*randn(num_steps, 1);
    ome_Q = sqrt(power_Q/on_SNR)*randn(num_steps, 1);
    varout(:, 1) = R + on_switch*ome_R;
    varout(:, 2) = S + on_switch*ome_S;
    varout(:, 3) = I + on_switch*ome_I;
    varout(:, 4) = Q + on_switch*ome_Q;
    varif = varout(num_steps-IF_length+1:num_steps, :);

% Multi-series causality analysis
    [t21, ~, ~, ~, ~] = multi_causality_est_all_new_tau(varif, 1, 1); 
    iscausal = abs(t21');
    iscausal(iscausal >= 0.1) = 1;
    iscausal(iscausal < 0.1) = 0;
    OmiR = OmiR + 1 - (iscausal(1, 2) + iscausal(1, 4) + iscausal(2, 3) + iscausal(2, 4) + iscausal(3, 2))/5;
    FasR = FasR + (iscausal(1, 3) + iscausal(2, 1) + iscausal(3, 1) + iscausal(4,1) + iscausal(4, 2) + iscausal(4, 3) + iscausal(3, 4))/5;
end
OmiR = OmiR/nsim;
FasR = FasR/nsim;
