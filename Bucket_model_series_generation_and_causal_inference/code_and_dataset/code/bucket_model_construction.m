function [R, S, I, Q] = bucket_model_construction(num_steps, pn_SNR, power_S, power_I, power_Q)
% 
% bucket_model_construction generates state variables R (rainfall), S (soil
% water content), I (interflow), and Q (runoff) using Eqn. (8) in this
% paper.
%
% Written by Yuhui Zhao (yuhuizhao73@foxmail.com)
%   
% On input:
%    num_steps: number of time steps for one simulation
%    pn_SNR: signal-noise-ratio of the process noise added to the system
%    power_S: power of the time series S
%    power_I: power of the time series I
%    power_Q: power of the time series Q
%
% On output:
%    R: rainfall
%    S: soil water content
%    I: interflow
%    Q: runoff
 
%% Specifications of the bucket model as those in Ombadi et al. (2020; Table 2)
S_max = 80; % maximum soil storage
S_0 = 40; % initial soil storage
Ks = 2.3; % speed of water depletion
delta = 10; % soil moisture amount below which no interflow occurs
xi = 0.6; % intensity of water depletion (between 0.5 and 1)
r = 0.8; % correlation coefficient of the process noise (red noise) between different time steps 
dt = 1; % time interval
 
%% Process noise calculation 
ome_S = sqrt(power_S/pn_SNR)*randn(num_steps,1);
ome_I = sqrt(power_I/pn_SNR)*randn(num_steps,1);
ome_Q = sqrt(power_Q/pn_SNR)*randn(num_steps,1);
ita_S = ome_S;
ita_I = ome_I;
ita_Q = ome_Q;
for t = 2:num_steps
    ita_S(t) = r*ita_S(t-1)+sqrt(1-r^2)*ome_S(t);
    ita_I(t) = r*ita_I(t-1)+sqrt(1-r^2)*ome_I(t);
    ita_Q(t) = r*ita_Q(t-1)+sqrt(1-r^2)*ome_Q(t);
end
 
%% Rainfall simulation
% Rainfall Occurence RO is simulated by a discrete-time first-order Markov chain model:
P_tilde = [0.8,0.2;0.2,0.8]; % transition probability [0-0,0-1;1-0,1-1]
mc = dtmc(P_tilde);
RO = simulate(mc,num_steps-1)-1;
% Rainfall amount Y is simulated by beta distribution:
B_alpha = 9;
B_beta = 6;
Y = betarnd(B_alpha,B_beta,num_steps,1)*80;
% Rainfall R = RO(t)*Y(t):
R = RO.*Y;
 
%% Construction of S, I,and Q
S = zeros(num_steps,1);
I = zeros(num_steps,1);
Q = zeros(num_steps,1);
S(1) = S_0;
for t = 2:num_steps
    S_delta = S(t-1)-delta;
    if S_delta > 0
       I(t) = Ks*S_delta^xi+ita_I(t);
    end
    S(t) = S(t-1)+dt*(R(t-1)-I(t-1)+ita_S(t));
    if S(t)<0
       S(t)=0;
    elseif S(t) >= S_max
       Q(t) = S(t-1)+R(t-1)-S_max+ita_Q(t);
       S(t) = S_max;
    end
end
end
