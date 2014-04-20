% This script serves as a wrapper code for implementing the Ensemble
% Kalman filter (EnKF). This script takes as inputs parameters those
% that were generated from the script: prob_model_wrapper.m, which
% is responsible for generating a synthetic true state and 
% corresponding measurements. If used with real data, then the 
% control parameter files should be generated appropriately and
% include the real measurements.  
%
% Author:   Bart Forman
% Created:  10 Feb 2012

clear all; tic

% Fix random number generation seed to that it is the same each execution
rng(1);

% Load necessary data generated from "true" simulation
load true_output.mat

% Number of pixels
N_pix=N_lat*N_lon

% Length of state vector
N_Y=N_y*N_pix

% Load nominal forcing inputs
load forcing.mat

% Specify number of replicates to use in ensemble
N_reps=50;

% Setup ICs (Y: states; Z:predicted meas.)
% Here a Gaussian distribution is used -- COULD USE OTHER DIST.!!
for k=1:N_reps
    y0(1:length(y0_bar),k)=mvnrnd(y0_bar,Cy0,1)';
end
Y(:,1,:) = y0;
N_Z = size(z_meas,1);
Z(1:N_Z,1,1:N_reps) = zeros(N_Z,1,N_reps);
index = 0;
jmeas = 1;

t=[];
% Loop through measurement intervals (i.e. between measurements)
for ii=1:length(z_meas)+1
    
    % Setup simulation control parameters
    if jmeas <= length(z_meas)
        if jmeas==1
            Day_beg=sim_Day_beg;
        else
            Day_beg=t_meas(jmeas-1);
        end
        Day_end=t_meas(jmeas);
    else
        Day_beg=t_meas(jmeas-1);
        Day_end=sim_Day_end;
    end

    disp(['Running measurement interval #: ' num2str(ii)])
    
    % Replicate Loop:
    for k=1:N_reps

        if k==1
            dims_Y=size(Y);
        end
        
        % Set control parameter
        control_par=[N_pix;N_Y;Day_beg;Day_end];

        % Set time-invariant parameters -- Here no uncertainty in
        % parameters is considered; if uncertainty is desired, could be
        % added here
        alpha=alpha;

        % Set IC: Either real IC at t=0 or IC after update
        y0_inp=Y(:,dims_Y(2),k);
        
        % Set forcing:
        % % Multiplicative error (constant in time)
        %forcing_inp = ...
        %    [forcing(:,1) ...
        %    forcing(:,2:end).*repmat(exp(mvnrnd(log_ubar_t,C_ut,1)),...
        %    length(forcing(:,1)),1)];
        % Multiplicative error (varying in time -- uncorrelated in time)
        forcing_inp = ...
            [forcing(:,1) ...
            forcing(:,2:end).*exp(mvnrnd(log_ubar_t,C_ut,length(forcing(:,1))))];
        
        % Call forward model -- PROPAGATION STEP
        [y_out,q_out,t_out] = ...
            forward_model(control_par,alpha,y0_inp,forcing_inp);
        
        if ii>1
            y_out=y_out(:,:,2:end);
        end
        
        % Compute predicted meas. at meas. time 
        z_out=H*reshape(y_out(:,:,end),N_Y,1);
        
        if k==1
            dims_y=size(y_out);
            dims_q=size(q_out);
            if ii>1
                dims_Q=size(Q);
            end
        end
        % Store state vector
        if ii==1
            Y(:,1:dims_y(3),k)=reshape(y_out,N_Y,dims_y(3));
            Q(:,1:dims_q(3),k)=reshape(q_out,3*N_pix,dims_q(3));
            Z(1:N_Z,ii,k)=z_out;
        else
            Y(:,dims_Y(2)+1:dims_Y(2)+dims_y(3),k)=reshape(y_out,N_Y,dims_y(3));
            Q(:,dims_Q(2)+1:dims_Q(2)+dims_q(3),k)=reshape(q_out,3*N_pix,dims_q(3));
            if ii <= length(z_meas)
                Z(1:N_Z,ii,k)=z_out;
            end
        end
        if k==1
            t=[t t_out(1:end)'];
        end
    end
    
    % Update replicates
    if jmeas <= length(z_meas)

        % Grabs states/predicted measurements at meas. time
        yrep([1:N_Y],[1:N_reps])=Y([1:N_Y],end,:);
        zrep([1:N_Z],[1:N_reps])=Z([1:N_Z],end,:);

        % Compute sample mean vectors
        ymean=mean(yrep,2)*ones(1,N_reps);
        zmean=mean(zrep)*ones(1,N_reps);

        % Compute sample covariance matrices
        Cyz=((yrep-ymean)*(zrep-zmean)')/(N_reps-1);
        Czz=((zrep-zmean)*(zrep-zmean)')/(N_reps-1);

        % compute Kalman gain
        kalm=Cyz/[Czz+Cvv];

        % Measurement error realizations
        v=mvnrnd(zeros(1,1),Cvv,N_reps)';

        % Update states
        index=size(Y,2);
        for k=1:N_reps
            yup(:,k)=yrep(:,k)+kalm*(z_meas(jmeas)+v(:,k)-zrep(:,k));
            Y(:,index+1,k)=yup(:,k);
            Z(:,index+1,k)=H*Y(:,index+1,k);
        end

        jmeas=jmeas+1;       

    end
end

% Reshape state vector 
dims=size(Y);
Y=reshape(Y,N_y,N_pix,dims(2),N_reps);
dims=size(Q);
Q=reshape(Q,3,N_pix,dims(2),N_reps);

save EnKF_output.mat t Y Q

toc
