% Script used to generate the "truth" for a single-state model of the
% form: y_k+1 = a * y_k + u_k
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

% Measurement model matrix (only one state included)
H = [1];

% Define initial state (only one state included)
y0 = [2];

% Define time-invariant parameters
alpha(1) = 0.9;

% TRUE Model forcing ("u" vector)
c1 = 0.04;
omega = 0.05;
for i=1:Nsteps
    utrue(i) = c1*sin(omega*i);
end
utrue = [0 utrue]; % add zero to the beginning of time series

% Call model (returns states over the entire simulation period)
[ytrue,ttrue] = Prob2_single_state_lin_model(y0,alpha,utrue,Nsteps,dt,t_beg);

% Define measurement vector (with normally distributed additive error)
j=1;
for i=5:5:Nsteps
    imeas(j) = i;
    tmeas(j) = ttrue(i);
    z(j) = H*ytrue(:,i) + randn*sqrt(Cvv);
    j = j+1;
end

% Save output for later use
save model_inp.mat t_beg Nsteps dt alpha utrue
save ytrue_single_state.mat ttrue ytrue
save zmeas_single_state.mat imeas tmeas z w_variance H Cvv
