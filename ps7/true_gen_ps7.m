% Script that specifies model inputs and calls model to generate true
% states and synthetic obs. for use in a DA scheme

clear all

% Fix random number generation seed to that it is the same each execution
rng(1);

% Number of time steps
time_fac = 1; % In case nonlinear model requires smaller time step
Nsteps = 400*time_fac;

% Time step
dt = 0.1/time_fac;

% simulation beginning time
t_beg = 0;

% Measurement error covariance
Cv = [0.05^2];

% Measurement model matrix (only measure state 2)
H = [1 0];

% Initial condition vector
y0true = [2; 0];

% Model parameters
f1 = 0.9;
f2 = 1;
f3 = 0;
f4 = 1;
c1 = 0.04;
omega = 0.05*150/Nsteps;

% Time-invariant parameter vector
alpha = [f1 f2 f3 f4 c1 omega];

% Model forcing
for i = 1:Nsteps
    u(:,i) = [c1*sin(omega*i);0];   % "b" vector
end

% Call model (returns states over simulation period)
[ytrue,ttrue] = lin_model(y0true,alpha,u,Nsteps,dt,t_beg);

% Define measurement vector (with normally distributed additive error)
j = 1;
for i = 30*time_fac:60*time_fac:Nsteps
    imeas(j) = i;
    tmeas(j) = ttrue(i);
    z(j) = H*ytrue(:,i) + randn*sqrt(Cv);
    j = j+1;
end

save model_inp.mat t_beg Nsteps dt alpha u
save ytrue.mat ttrue ytrue
save zmeas.mat imeas tmeas z Cv H

%% Plot output
figure(1)
clf
subplot(2,1,1)
plot(ttrue,ytrue(1,:),'k','LineWidth',2)
hold on
ylabel('y_1','FontSize',12);xlabel('t','FontSize',12)
grid

subplot(2,1,2)
plot(ttrue,ytrue(2,:),'k','LineWidth',2)
hold on
plot(tmeas,z,'mo','LineWidth',2)
ylabel('y_2','FontSize',12);xlabel('t','FontSize',12)
grid
legend('true','meas.')
