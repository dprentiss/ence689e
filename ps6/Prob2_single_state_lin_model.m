% Script used to run a single-state linear model for use
% in a Kalman Filter (KF) application.
%
% Author:   Bart Forman

function [y,t] = Prob2_single_state_lin_model(y0,alpha,u,Nsteps,dt,t_beg)

% Set IC
y(:,1) = y0;

% Set parameters
A = alpha(1);

% Set forcing
b = dt*u;

% Run forward model. NOTE: time indexing starts at "1" since MATLAB
% does not support "0" indices
t(1)=t_beg;
for i=1:Nsteps
    
    t(i+1) = t_beg+i*dt; % time vector

    % propagation of state to next time step
    y(:,i+1) = A*y(:,i) + b(:,i);             

end
