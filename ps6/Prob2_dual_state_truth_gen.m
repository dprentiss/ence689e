% Script used to generate the "truth" for a two-state system for use
% in a Kalman Filter (KF) application.
%
% Author:   Bart Forman

clear all; clc;

% Fix random number generation seed to that it is the same each execution
rng(1);

% Number of time steps
Nsteps=150;

% Time step
dt=1;

% Simulation beginning time
t_beg=0;

% White noise variance for w
w_variance = 0.10^2;

% Measurement error (co)variance
Cvv = 0.05^2;

% Measurement model matrix (only measures state 1)
H = [1 0];

% Define initial state
y0true = [2 0];

% Define time-invariant parameter vector
alpha(1) = 0.9;
alpha(2) = 1.0;
alpha(3) = 0;
alpha(4) = 0.50; % Q: What is the "best" parameter to use here? (HINT: It is not 0.5)
alpha(5) = 1;
alpha(6) = 0;
alpha(7) = 0;
alpha(8) = 1;
c1 = 0.04;
omega = 0.05;
alpha(9) = c1;
alpha(10) = omega;

% TRUE Model forcing ("u" vector)
for i=1:Nsteps
    utrue(i)=c1*sin(omega*i);
end
utrue = [0 utrue]; % add zero to the beginning of time series

% Call model (returns states over the entire simulation period including the
% estimated forcing)
[ytrue,ttrue] = ...
    Prob2_dual_state_lin_model(y0true,alpha,utrue,Nsteps,dt,t_beg,w_variance);

% Add "true" forcing to second state in the "true" augmented state vector
ytrue(2,:) = utrue;

% Define measurement vector (with normally distributed additive error)
j=1;
for i=5:10:Nsteps
    imeas(j) = i;
    tmeas(j) = ttrue(i);
    z(j) = H*ytrue(:,i) + randn*sqrt(Cvv);
    j = j+1;
end

% Save output for later use (NOTE: measurements were generate in the single-
% state simulation and are not repeated here for simplicity)
save model_inp.mat t_beg Nsteps dt alpha utrue
save ytrue.mat ttrue ytrue
save zmeas.mat imeas tmeas z w_variance H Cvv
