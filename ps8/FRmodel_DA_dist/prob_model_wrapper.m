% This script serves as a wrapper for calling a forward model for use in 
% probabilistic open-loop modeling and/or synthetic data assimilation
% studies. It is meant to be modular (i.e., it could be used with other
% models. Here it is applied to the Force-Restore land-surface model.
%
% Author:   Bart Forman
% Created:  10 Feb 2012

clear all; tic

% Fix random number generation seed to that it is the same each execution
rng(1);

% Basic model inputs:
N_y = 4;  % # of states per pixel (y = [Ts, Td, W1, W2]^T at each pixel)

% Size of domain
N_lat = 1     % number or pixels in y-direction
N_lon = 1     % number of pixels in x-direction

% Total number of pixels
N_pix = N_lat*N_lon

% Length of state vector
N_Y = N_y*N_pix

% Simulation control parameters
sim_Day_beg = 160.;
sim_Day_end = 161.;

% Specify nominal (mean) inputs and DEM information. NOTE: The DEM is 
% 151x151, hence the domain can be up to this size. If a smaller
% region is specified, the wrapper uses the upper-left corner.
load DEM_data.mat   
elev_inp = reshape(elev(1:N_lat,1:N_lon),N_pix,1);
slope_inp = reshape(slope(1:N_lat,1:N_lon),N_pix,1);
aspect_inp = reshape(aspect(1:N_lat,1:N_lon),N_pix,1);
lat_inp = reshape(repmat(lat_deg(1:N_lat)',[1 N_lon]),N_pix,1);

% Load nominal forcing inputs
load forcing.mat
% Specify forcing uncertainty model
% Here a multiplicative lognormal error model is used
cov_ut = [0 0 0 0 0.25 0 0];    % Coeff. of variation for each forcing var.
                                % order corresponds to forcing vector
                                % elements; here only precip. uncertainty
                                % considered
% Transform to parameters needed for mvnrnd -- Note: this assumes the
% forcing variables are uncorrelated with each other
var_u = log(1+cov_ut.^2);         % Compute variance     
C_ut = diag(var_u);               % Compute covariance
log_ubar_t = log(1)-0.5*var_u;    % Mean parameter for mvnrnd

% Set initial conditions
% Surface/deep temp.:
Ts0 = 294;
Td0 = Ts0+1;

% Surface/deep soil moisture:
W10 =  0.5;
W20 = 0.5;
y0_bar  =  ...
    [Ts0*ones(1,N_pix); ...
    Td0*ones(1,N_pix); ...
    W10*ones(1,N_pix); ...
    W20*ones(1,N_pix)];

% Mean vector
y0_bar = reshape(y0_bar,N_Y,1);

% Assume states ICs are uncorrelated with each other (could change)
% Construct covariance matrix
% Covariance between states at a single pixel
Cy0_pix = [1^2 0   0      0; ...
           0   1^2 0      0; ...
           0   0   0.05^2 0; ...
           0   0   0      0.05^2];
% Construct for full state vector (i.e. over spatial domain)
Cy0 = [];
for i = 1:N_pix
    Cy0 = blkdiag(Cy0,Cy0_pix);
end

% Get some info. from user
run_flag = ...
    input(['Enter "0" for "true" realization simulation;',...
    ' "1" for probabilistic open-loop simulation:  ']);
if (run_flag==0)
    N_reps = 1;
    % Define measurements the average surface soil moisture over
    % the entire domain.
    
    % Measurement error covariance
    Cvv = [0.03^2];
    
    % Measurement model matrix for a single pixel
    % NOTE: This does NOT need to be linear, but it is applied so here.
    H_pix = [0 0 1 0]; % Measure surface soil moisture only
    
    % Note: This takes the average over the N_pix pixels (i.e. from a
    % coarse satellite measurement)
    H = repmat(H_pix,1,N_pix)/N_pix;
    
    % Define meas. times in a given day (day fraction):
    daily_meas_times = [10*60*60/86400];    % 10 am
    t_meas = [];
    for i_day = sim_Day_beg:sim_Day_end-1
        t_meas_day = i_day+daily_meas_times;
        t_meas = [t_meas t_meas_day'];
    end
else
    % Number of replicates (ensemble size)
    N_reps = input('Enter number of replicates:  ');
end

for i_rep = 1:N_reps
    
    disp(['Running replicate #: ' num2str(i_rep)])

    % Set control parameter
    control_par = [N_pix; N_Y; sim_Day_beg; sim_Day_end];
    
    % Sample from IC pdfs
    if run_flag==0 % Make the truth different from the open-loop mean
        y0_bar = ...
            reshape([296*ones(1,N_pix); 297*ones(1,N_pix); ...
            0.4*ones(1,N_pix); 0.4*ones(1,N_pix)],N_Y,1);
        y0 = mvnrnd(y0_bar,Cy0)';
    else
        y0 = mvnrnd(y0_bar,Cy0)';
    end
    
    % Set time-invariant parameters
    alpha = [elev_inp slope_inp aspect_inp lat_inp];

    % Apply multiplicative, lognormal perturbations (constant in time)
    if run_flag==0
        forcing_inp = ...
            [forcing(:,1) ...
            forcing(:,2:end).*repmat([1 1 1 1 0.5 1 1], ...
            length(forcing(:,1)),1)]; % hardwire low-bias precip.
    else
        forcing_inp = ...
            [forcing(:,1) ...
            forcing(:,2:end).*repmat(exp(mvnrnd(log_ubar_t,C_ut,1)),...
            length(forcing(:,1)),1)];
    % % Option below would add a different multiplicative error to each time
    %forcing_inp = ...
    %    [forcing(:,1) ...
    %    forcing(:,2:end).*exp(mvnrnd(log_ubar_t,C_ut,...
    %    length(forcing(:,1))))];
    end
    
    % Call forward model
    [y_out,q_out,t_out] = forward_model(control_par,alpha,y0,forcing_inp);

    % Save states and fluxes for analysis
    Y(:,:,:,i_rep) = y_out;
    Q(:,:,:,i_rep) = q_out;

end


if (run_flag==0)
    % Save true states/measurements and control parameters, etc.
    dims = size(Y);
    YY = reshape(Y,N_Y,dims(3));
    z_true = interp1(t_out,(H*YY),t_meas);
    z_meas = z_true + randn(length(t_meas),1)'*sqrt(Cvv);
    t_true = t_out;
    for nn = 1:length(t_meas)
        imeas(nn) = find(abs(t_true-t_meas(nn))<10^-5);
    end
    Y_true = Y;
    Q_true = Q;
    save true_output.mat N_y N_lat N_lon sim_Day_beg sim_Day_end alpha log_ubar_t C_ut y0_bar Cy0 t_true Y_true Q_true H z_true z_meas Cvv t_meas

else
    % Save open-loop outputs
    t_open = t_out;
    Y_open = Y;
    Q_open = Q;    
    save prob_open_loop_output.mat t_open Y_open Q_open
end

toc
