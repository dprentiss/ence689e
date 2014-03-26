% This script applies the Kalman Filter (KF) to solve the state 
% estimation problem using the same two-state linear model used
% with the linear model. The measurements and error statistics
% are loaded and subsequently used in the algorithm.
%
% Author:   Bart Forman

clear all; clc; tic

% Load inputs
load model_inp.mat
load ytrue.mat
load zmeas.mat

% Define the initial state vector
y0bar = [ 2; 0 ];

% Specify prior mean and initial error covariance
ybar = y0bar;
Cy0y0 = [0 0; 0 0]; % Q: What is the initial error covariance (HINT: It is not all zeros)
Cyy = Cy0y0;
t = t_beg;

% Specify white noise covariance
Cww = [0 0; 0 w_variance];

%%% Kalman Filter to Estimate Posterior %%%%

jmeas = 1;
    
% Loop through measurement intervals
for ii = 1:length(z)+1
    if jmeas <= length(z)
        meastime = imeas(jmeas);
    end

    if jmeas==1
        ind1 = 1;
    else
        ind1 = ind2+1+1;
    end
    if jmeas <= length(z) 
        ind2 = imeas(jmeas)-1 +(jmeas-1);
    else
        ind2 = Nsteps+length(z);
    end

    % PROPAGATION STEP:
    % Determine proper nsteps
    nsteps = ind2-ind1+1;
    % Grab proper u
    if jmeas==1 
        uind1 = 1; uind2 = imeas(jmeas)-1;
    elseif jmeas <=length(z)
        uind1 = imeas(jmeas-1); uind2 = imeas(jmeas)-1;
    else
        uind1 = imeas(jmeas-1); uind2 = Nsteps;
    end
    uu = utrue(:,uind1:uind2);

    % Assign proper initial mean vector/covariance
    if jmeas>1
        t_beg = t(end);
        y0bar = ybar(:,end);
        Cy0y0 = Cyy(:,:,end);
    end
    [ybar_out,Cy_out,t_out] = ...
        Prob2_open_loop_lin_model(y0bar,Cy0y0,alpha,uu,nsteps,dt,t_beg,Cww);

    % Assign output variables
    t = [t t_out(2:end)];
    ybar = [ybar ybar_out(:,2:end)];
    j = size(Cyy,3);P = size(Cy_out,3);
    for p=2:P
        Cyy(:,:,j+1) = Cy_out(:,:,p);    
        j = j+1;
    end

    % Plot output
    figure(2)
    clf
    subplot(2,1,1)
    plot(ttrue,ytrue(1,:),'k','LineWidth',2)
    hold on
    plot(tmeas,z,'mo','LineWidth',2)
    plot(t,ybar(1,:),'b','LineWidth',2)
    plot(t,squeeze(ybar(1,:))'+2*squeeze(sqrt(Cyy(1,1,:))),'b:','LineWidth',2)
    plot(t,squeeze(ybar(1,:))'-2*squeeze(sqrt(Cyy(1,1,:))),'b:','LineWidth',2)
    xlabel('t','FontWeight','bold'); ylabel('y_1','FontWeight','bold');grid
    legend('Truth','Measurements','Posterior','+/- 2 std. dev.','Location','NorthEast')
    title('Propagation Step','FontWeight','bold')

    subplot(2,1,2)
    plot(ttrue,ytrue(2,:),'k','LineWidth',2)
    hold on
    plot(t,ybar(2,:),'b','LineWidth',2)
    plot(t,squeeze(ybar(2,:))'+2*squeeze(sqrt(Cyy(2,2,:))),'b:','LineWidth',2)
    plot(t,squeeze(ybar(2,:))'-2*squeeze(sqrt(Cyy(2,2,:))),'b:','LineWidth',2)
    xlabel('t','FontWeight','bold'); ylabel('y_2','FontWeight','bold');grid

    % UPDATE STEP
    if jmeas <= length(z)
        Cyz = Cyy(:,:,end)*H';
        Czz = H*Cyy(:,:,end)*H'+Cvv;
        K = Cyz*inv(Czz);
        index = size(ybar,2);
        ybar(:,index+1) = ybar(:,index) + K * (z(:,jmeas) - H*ybar(:,index)) ;
        Cyy(:,:,index+1) = [eye(2) - K*H]*Cyy(:,:,index);
        t = [t t(end)];     
        jmeas = jmeas+1;
    end

    if jmeas<=length(z)
        % Plot output
        figure(2)
        clf
        subplot(2,1,1)
        plot(ttrue,ytrue(1,:),'k','LineWidth',2)
        hold on
        plot(tmeas,z,'mo','LineWidth',2)
        plot(t,ybar(1,:),'b','LineWidth',2)
        plot(t,squeeze(ybar(1,:))'+2*squeeze(sqrt(Cyy(1,1,:))),'b:','LineWidth',2)
        plot(t,squeeze(ybar(1,:))'-2*squeeze(sqrt(Cyy(1,1,:))),'b:','LineWidth',2)
        xlabel('t','FontWeight','bold'); ylabel('y_1','FontWeight','bold'); grid
        title('Update Step','FontWeight','bold')
        legend('Truth','Measurements','Posterior','+/- 2 std. dev.','Location','NorthEast')
        subplot(2,1,2)
        plot(ttrue,ytrue(2,:),'k','LineWidth',2)
        hold on
        plot(t,ybar(2,:),'b','LineWidth',2)
        plot(t,squeeze(ybar(2,:))'+2*squeeze(sqrt(Cyy(2,2,:))),'b:','LineWidth',2)
        plot(t,squeeze(ybar(2,:))'-2*squeeze(sqrt(Cyy(2,2,:))),'b:','LineWidth',2)
        xlabel('t','FontWeight','bold'); ylabel('y_2','FontWeight','bold'); grid        
    end

end
    
toc

%% Plot the KF output at the final time step

%%%%% add code here (HINT: Reuse the code above accordingly) %%%%%

%%%%% Also, add the following title above the first state estimate: %%%%%
%title('Final Estimate','FontWeight','bold')
