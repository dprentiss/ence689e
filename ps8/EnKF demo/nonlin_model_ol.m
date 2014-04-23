% nonlinear model for state propagation

function [Y,t] = nonlin_model(Y0,alpha,u,Nsteps,dt,t_beg)

% Set IC
Y(:,1) = Y0;

% Set parameters
f1 = alpha(1); f2 = alpha(2); f3 = alpha(3); f4 = alpha(4);
%A = [1+f1*dt f2*dt; f3*dt 1+f4*dt];   % "A" matrix 
%c1 = alpha(5); omega = alpha(6);

% Set forcing
b = dt*u;

% Run forward model. Note that time indexing starts at "1" 
% since MATLAB does not support "0" indices
t(1) = t_beg;
for i = 1:Nsteps
    
    t(i+1) = t_beg+i*dt; % time vector
    AA = [1+f1*dt Y(3,i)*dt 0; f3*dt 1+f4*dt 0; 0 0 1];
    % propagation of state to next time step
    Y(:,i+1) =  AA*Y(:,i) + [b(:,i);0];             

end

