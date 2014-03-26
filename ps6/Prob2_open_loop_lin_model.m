function [ybar,Cyy,t] = Prob2_open_loop_lin_model(y0bar,Cy0y0,alpha,u,Nsteps,dt,t_beg,Cww)

% Set initial conditions
ybar(:,1)=y0bar;
Cyy(:,:,1)=Cy0y0;

% Establish matrix A
A = [alpha(1) alpha(2); alpha(3) alpha(4)];

% Establish matrix C
c = [alpha(5) alpha(6); alpha(7) alpha(8)];

% Set AR(1) forcing model error
for i=1:Nsteps
    w(1,i) = 0;
    w(2,i) = mvnrnd(0,Cww(2,2),1);
end

% Run forward model. NOTE: time indexing starts at "1" since MATLAB
% does not support "0" indices
t(1)=t_beg;
for i=1:Nsteps
    
    t(i+1)=t_beg+i*dt; % time vector
    
    % Set forcing
    ybar(2,i) = u(1,i);
    
    % Propagate mean state vector to next time step
    ybar(:,i+1) = A*ybar(:,i);
    
    % Propagate error covariance matrix to next time step
    Cyy(:,:,i+1) = A*Cyy(:,:,i)*A' + c*Cww*c';
    
end
