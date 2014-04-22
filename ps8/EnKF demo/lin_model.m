% linear model for state propagation

function [y,t] = lin_model(y0,alpha,u,Nsteps,dt,t_beg)

% Set IC
y(:,1) = y0;

% Set parameters
f1 = alpha(1); f2 = alpha(2); f3 = alpha(3); f4 = alpha(4);
A = [1+f1*dt f2*dt; f3*dt 1+f4*dt];   % "A" matrix 
%c1 = alpha(5); omega = alpha(6);

% Set forcing
b = dt*u;

% Run forward model. Note that time indexing starts at "1" 
% since MATLAB does not support "0" indices
t(1) = t_beg;
for i = 1:Nsteps
    
    t(i+1) = t_beg+i*dt; % time vector

    % propagation of state to next time step
    y(:,i+1) =  A*y(:,i) + b(:,i);             

end
