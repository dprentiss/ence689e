% Script used to run a linear model for use in a Kalman Filter (KF)
% application.
%
% Author:   Bart Forman

function [y,t] = Prob2_dual_state_lin_model(y0,alpha,u,Nsteps,dt,t_beg,w_variance)

% Set IC
y(:,1) = y0;

% Establish matrix A
A = [alpha(1) alpha(2); alpha(3) alpha(4)];

% Establish matrix C
c = [alpha(5) alpha(6); alpha(7) alpha(8)];

% Set AR(1) forcing model error
for i=1:Nsteps
    w(1,i) = 0;
    w(2,i) = mvnrnd(0,w_variance,1);
end

% Run forward model. NOTE: time indexing starts at "1" since MATLAB
% does not support "0" indices
t(1)=t_beg;
for i=1:Nsteps
    
    t(i+1) = t_beg+i*dt; % time vector
    
    % Set forcing
    y(2,i) = u(1,i);
    
    % Propagate mean state vector to next time step
    y(:,i+1) = A*y(:,i) + c*w(:,i);
    
end
