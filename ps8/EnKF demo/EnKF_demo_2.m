% EnKF.m
%
% This script applys the Ensemble Kalman Filter (EnKF) to solve
% the augmented state estimation problem using the three-state 
% nonlinear model. The measurements and error statistics
% are loaded and used in the algorithm.

clear all; tic

% Load inputs
load model_inp.mat
load ytrue.mat
load zmeas.mat

% Fix random number generation seed to that it is the same each execution
rng(1);

%f2 = -0.1875;
%alpha(2) = f2;

% Augmented measurement model
H = [0 1 0];

%% Specify prior mean/covariance at inital time
ybar = [0;0];
Cy0y0 = [4 1.5; 1.5 1];

f2bar = -0.1875;
Cf2 = 0.2^2;

Y0bar = [ybar;f2bar];
Cy0y0 = [Cy0y0(1,:) 0; Cy0y0(2,:) 0; 0 0 Cf2];
Ybar = [ybar;f2bar];
CYY = Cy0y0;

% Specify number of replicates to use (ensemble size)
Nreps = 100;

%% Run EnKF
% Setup ICs (y: states; zz:predicted meas.)
for k = 1:Nreps
    Y0(1:length(Ybar),k) = mvnrnd(Ybar,Cy0y0,1)';
    Y0(3,k) = -abs(Y0(3,k));
    zz(1,1,k) = H*Y0(:,k);
end
Y(:,1,:) = Y0;
zz(:,1,:) = H*Y0;
t = [t_beg];
index = 0;
jmeas = 1;
    
% Loop through measurement intervals (i.e. between measurements)
for ii = 1:length(z)+1
    if jmeas <=  length(z)
        meastime = imeas(jmeas);
    end

    % Loop through replicates
    for k = 1:Nreps

        % Define appropriate indices for running model between meas.
        % ind1: starting index; ind2:ending index
        if k == 1
            if jmeas == 1
                ind1 = 1;
            else
                ind1 = ind2+1+1;
            end
            if jmeas <=  length(z)
                ind2 = imeas(jmeas)-1 +(jmeas-1);
            else
                ind2 = Nsteps+length(z);
            end
        end

        % Call model
        % Determine proper nsteps
        nsteps = ind2-ind1+1;
        % Grab proper forcing vector u
        % uind1: starting index for u; uind2: ending index for u
        if jmeas == 1 
            uind1 = 1; uind2 = imeas(jmeas)-1;
        elseif jmeas <= length(z)
            uind1 = imeas(jmeas-1); uind2 = imeas(jmeas)-1;
        else
            uind1 = imeas(jmeas-1); uind2 = Nsteps;
        end
        uu = u(:,uind1:uind2);

        % Assign proper initial conditions 
            t_beg = t(end);
            Y0_inp = Y(:,index+1,k);

        % Call model (returns states over simulation period)
        [Y_out,t_out] = nonlin_model_2(Y0_inp,alpha,uu,nsteps,dt,t_beg);

        % Compute predicted measurements
        z_out = H*Y_out;

        % Put model outputs/meas. into appropriate vectors
        filly = [Y_out(:,2:end)];
        fillz = [z_out(:,2:end)];
        if k == 1
            ind1y = size(Y(:,:,k),2);
        end
        ind2y = size(filly,2);
        Y(:,ind1y+1:ind1y+ind2y,k) = filly;
        zz(:,ind1y+1:ind1y+ind2y,k) = fillz;
        if k == 1
            t = [t t_out(2:end)];
        end

    end
    figure(5)
    clf
    subplot(3,1,1)
    plot(ttrue,ytrue(1,:),'k','LineWidth',2)
    hold on
    plot(t,squeeze(mean(Y(1,:,:),3)),'b','LineWidth',2)
    plot(t,squeeze(Y(1,:,:)),'c:','LineWidth',2)
    plot(ttrue,ytrue(1,:),'k','LineWidth',2)
    plot(t,squeeze(mean(Y(1,:,:),3)),'b','LineWidth',2)
    xlabel('t','FontSize',12); ylabel('y_1','FontSize',12);grid
    title('Propagation Step','FontWeight','bold')

    subplot(3,1,2)
    plot(ttrue,ytrue(2,:),'k','LineWidth',2)
    hold on
    plot(tmeas,z,'mo','LineWidth',2)
    plot(t,squeeze(mean(Y(2,:,:),3)),'b','LineWidth',2)
    plot(t,squeeze(Y(2,:,:)),'c:','LineWidth',2)
    legend('Truth','Measurements','Ensemble Mean','Replicates')
    plot(ttrue,ytrue(2,:),'k','LineWidth',2)
    plot(t,squeeze(mean(Y(2,:,:),3)),'b','LineWidth',2)
    plot(tmeas,z,'mo','LineWidth',2)
    xlabel('t','FontSize',12); ylabel('y_2','FontSize',12);grid

    subplot(3,1,3)
    plot([0 Nsteps*dt], [-2*0.1875 -2*0.1875],'k-','LineWidth',2)
    hold on
    plot(t,squeeze(mean(Y(3,:,:),3)),'b','LineWidth',2);hold on
    plot(t,squeeze(Y(3,:,:)),'c:','LineWidth',2)
    plot([0 Nsteps*dt], [-2*0.1875 -2*0.1875],'k-','LineWidth',2)
    plot(t,squeeze(mean(Y(3,:,:),3)),'b','LineWidth',2);hold on
    axis([0 Nsteps*dt -1 1])
    xlabel('t','FontSize',12); ylabel('f_2','FontSize',12);grid
    
    pause
    
    % Update replicates
    if jmeas <=  length(z)
        % Grabs states/predicted measurements at meas. time
        ylast([1:3],[1:Nreps]) = Y([1:3],end,:);
        yrep = ylast;
        zlast([1:1],[1:Nreps]) = zz([1:1],end,:);
        zrep = zlast;
        % Compute sample mean vectors
        ymean = mean(yrep,2)*ones(1,Nreps);
        zmean = mean(zrep)*ones(1,Nreps);
        % Compute sample covariance matrices
        Cyz = ((yrep-ymean)*(zrep-zmean)')/(Nreps-1);
        Czz = ((zrep-zmean)*(zrep-zmean)')/(Nreps-1);
        % Compute Kalman gain
        kalm = Cyz/[Czz+Cvv];
        % Measurement error realizations
        v = mvnrnd(zeros(1,1),Cvv,Nreps)';
        % Update states
        index = size(Y,2);
        for k = 1:Nreps
            yup(:,k) = yrep(:,k)+kalm*(z(jmeas)+v(:,k)-zrep(:,k));
            Y(:,index+1,k) = yup(:,k);
            Y(3,index+1,k) = -abs(Y(3,index+1,k));
            zz(:,index+1,k) = H*Y(:,index+1,k);
        end
        % Add extra time element at measurement (one before; one after
        % update)
        t = [t t(end)];

        jmeas = jmeas+1;

    end

    if jmeas <= length(z)
    figure(5)
    clf
    subplot(3,1,1)
    plot(ttrue,ytrue(1,:),'k','LineWidth',2)
    hold on
    plot(t,squeeze(mean(Y(1,:,:),3)),'b','LineWidth',2)
    plot(t,squeeze(Y(1,:,:)),'c:','LineWidth',2)
    plot(ttrue,ytrue(1,:),'k','LineWidth',2)
    plot(t,squeeze(mean(Y(1,:,:),3)),'b','LineWidth',2)
    xlabel('t','FontSize',12); ylabel('y_1','FontSize',12);grid
    title('Update Step','FontWeight','bold')

    subplot(3,1,2)
    plot(ttrue,ytrue(2,:),'k','LineWidth',2)
    hold on
    plot(tmeas,z,'mo','LineWidth',2)
    plot(t,squeeze(mean(Y(2,:,:),3)),'b','LineWidth',2)
    plot(t,squeeze(Y(2,:,:)),'c:','LineWidth',2)
    legend('Truth','Measurements','Ensemble Mean','Replicates')
    plot(ttrue,ytrue(2,:),'k','LineWidth',2)
    plot(t,squeeze(mean(Y(2,:,:),3)),'b','LineWidth',2)
    plot(tmeas,z,'mo','LineWidth',2)
    xlabel('t','FontSize',12); ylabel('y_2','FontSize',12);grid

    subplot(3,1,3)
    plot([0 Nsteps*dt], [-2*0.1875 -2*0.1875],'k-','LineWidth',2)
    hold on
    plot(t,squeeze(mean(Y(3,:,:),3)),'b','LineWidth',2);hold on
    plot(t,squeeze(Y(3,:,:)),'c:','LineWidth',2)
    plot([0 Nsteps*dt], [-2*0.1875 -2*0.1875],'k-','LineWidth',2)
    plot(t,squeeze(mean(Y(3,:,:),3)),'b','LineWidth',2);hold on
    axis([0 Nsteps*dt -1 1])
    xlabel('t','FontSize',12); ylabel('f_2','FontSize',12);grid

    %pause

    end

end

toc

% Plot output
figure(5)
clf
subplot(3,1,1)
plot(ttrue,ytrue(1,:),'k','LineWidth',2)
hold on
plot(t,squeeze(mean(Y(1,:,:),3)),'b','LineWidth',2)
plot(t,squeeze(Y(1,:,:)),'c:','LineWidth',2)
plot(ttrue,ytrue(1,:),'k','LineWidth',2)
plot(t,squeeze(mean(Y(1,:,:),3)),'b','LineWidth',2)
xlabel('t','FontSize',12); ylabel('y_1','FontSize',12);grid
title('Final Estimate','FontWeight','bold')

subplot(3,1,2)
plot(ttrue,ytrue(2,:),'k','LineWidth',2)
hold on
plot(tmeas,z,'mo','LineWidth',2)
plot(t,squeeze(mean(Y(2,:,:),3)),'b','LineWidth',2)
plot(t,squeeze(Y(2,:,:)),'c:','LineWidth',2)
legend('Truth','Measurements','Ensemble Mean','Replicates')
plot(ttrue,ytrue(2,:),'k','LineWidth',2)
plot(t,squeeze(mean(Y(2,:,:),3)),'b','LineWidth',2)
plot(tmeas,z,'mo','LineWidth',2)
xlabel('t','FontSize',12); ylabel('y_2','FontSize',12);grid

subplot(3,1,3)
plot([0 Nsteps*dt], [-2*0.1875 -2*0.1875],'k-','LineWidth',2)
hold on
plot(t,squeeze(mean(Y(3,:,:),3)),'b','LineWidth',2);hold on
plot(t,squeeze(Y(3,:,:)),'c:','LineWidth',2)
plot([0 Nsteps*dt], [-2*0.1875 -2*0.1875],'k-','LineWidth',2)
plot(t,squeeze(mean(Y(3,:,:),3)),'b','LineWidth',2);hold on
axis([0 Nsteps*dt -1 1])
xlabel('t','FontSize',12); ylabel('f_2','FontSize',12);grid
        