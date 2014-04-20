% This is an intermediate script designed to run a point-scale land surface
% model in a spatially distributed mode.  The control parameters,
% time-invariant model parameters, initial conditions, and forcing are
% passed to this routine, which then saves for each pixel in the domain and
% then calls the single-pixel force-restore land surface model.  The states
% and fluxes are returned to the calling script

function [y,q,t]=forward_model(control_par,alpha,y0,forcing)

N_pix=control_par(1);
N_Y=control_par(2);
Day_beg=control_par(3);
Day_end=control_par(4);

elev_inp=alpha(:,1);
slope_inp=alpha(:,2);
aspect_inp=alpha(:,3);
lat_inp=alpha(:,4);

elev_nom=min(elev_inp);   % meters above sea-level

y0_mat=reshape(y0,N_Y/N_pix,N_pix);

% Run a spatially distributed simulation (this assumes the model is a
% single pixel model that is not already vectorized -- model called/run
% sequentially

for i=1:N_pix

    %disp(['Running pixel #: ' num2str(i)])
    
    % For simplicity assume all time-invariant parameters are uniform in
    % space and only forcing is different
    
    % These next three function calls are simply writing the necessary
    % input files for the model
    
    % 1)  Write out control parameter file [can pass arguments that vary]
    Ts0=y0_mat(1,i) ; Td0=y0_mat(2,i) ; W10=y0_mat(3,i); W20=y0_mat(4,i);
    save_control_params(Ts0,Td0,W10,W20,Day_beg,Day_end);
    
    % 2) Write out time-invariant parameter file [could pass arguments that vary]
    save_time_invar_params(elev_inp(i),slope_inp(i),aspect_inp(i),lat_inp(i),elev_nom);
    
    % 3) Write out forcing file [could pass forcing that varies]
    save_forcing_file(forcing);
    
    % Call single-pixel Force-Restore Model
    FRmodel
    
    % Load output file and collect state/flux vector
    load states.out
    load fluxes.out
    
    % Store states and fluxes
    t=states(:,1);
    y(:,i,:)=states(:,2:end)';
    q(:,i,:)=fluxes(:,2:end)';
    
end
