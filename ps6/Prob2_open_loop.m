% This script runs a dual state, linear model without assimilation.
%
% Author:   Bart Forman

clear all; clc; tic

% Load inputs
load model_inp.mat
load ytrue.mat
load zmeas.mat

% Specify initial conditions for the 2-state system
y0bar = [ 2; 0 ];
ybar_FULL = y0bar;

% Specify error covariance for the 2-state system at inital time
Cy0y0 = [0 0; 0 0]; % Q: What is the initial error covariance (HINT: It is not all zeros)

% Initialize the mean vector and covariance matrix
ybar = y0bar;
Cyy = Cy0y0;
t = t_beg;

% Define white noise sequence
w = sqrt(w_variance)*randn([Nsteps 1]);

% Specify the error covariance matrix for the forcing
Cww = [0 0; 0 w_variance];

% Define forcing based on AR(1) model
u_AR(1) = utrue(1); % set AR forcing equal to utrue forcing in the beginning
rho = alpha(4);
for i=1:1:Nsteps
    u_AR(i+1) = rho * utrue(i) + w(i);
end

%% Open-loop Simulation
for t = 1:1:Nsteps
    
    % Store time step information for later use
    t_beg = t(end);
    
    % Assign proper initial mean vector/covariance
    if t == 1
        y0bar = ybar(:,end);
        Cy0y0 = Cyy(:,:,end);
    else
        y0bar = ybar_out(:,end);
        Cy0y0 = Cy_out(:,:,end);
    end
    
    % Define number of steps to propagate the model forward in time
    nsteps = 1;
    
    % Propagate the mean and covariance forward in time
    [ybar_out,Cy_out,t_out] = ...
        Prob2_open_loop_lin_model(y0bar,Cy0y0,alpha,u_AR(t),nsteps,dt,t_beg,Cww);
    
    % Capture output variables
    ybar_FULL(:,t+1) = ybar_out(:,end);
    Cyy_FULL(:,:,t+1) = Cy_out(:,:,end);
    
end

