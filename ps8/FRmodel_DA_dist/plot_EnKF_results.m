% Script use to plot open-loop results relative to the "truth".
%
% Author:   Bart Forman
% Created:  10 Feb 2012

clear all; tic;

load true_output.mat
load prob_open_loop_output.mat
load EnKF_output.mat

[N,M,P]=size(Y_true);

% Show some single pixel results
pix_num=input(['Enter pixel number for time-series plotting (1...' num2str(M) '): ']);

for i_state=1:4
    figure(i_state+1); clf
    
    subplot(1,2,1)
    hold on
    h1 = plot(t_open,squeeze(Y_open(i_state,pix_num,:,:)),':c','LineWidth',2);
    h2 = plot(t_true,squeeze(Y_true(i_state,pix_num,:)),'-k','LineWidth',2);
    h3 = plot(t_open,squeeze(mean(Y_open(i_state,pix_num,:,:),4)),'-m','LineWidth',2);
    xlabel('DOY')
    ylabel(['Y_' num2str(i_state)])
    legend([h2; h3; h1(1)],'Truth','OL Mean','OL Replicates')
    grid on
    
    subplot(1,2,2)
    hold on
    h1 = plot(t,squeeze(Y(i_state,pix_num,:,:)),':c','LineWidth',2);
    h2 = plot(t_true,squeeze(Y_true(i_state,pix_num,:)),'-k','LineWidth',2);
    h3 = plot(t,squeeze(mean(Y(i_state,pix_num,:,:),4)),'-m','LineWidth',2);
    xlabel('DOY')
    legend([h2; h3; h1(1)],'Truth','EnKF Mean','EnKF Replicates')
    grid on
end

for i_state=1:3
    figure(i_state+5); clf
    
    subplot(1,2,1)
    hold on
    h1 = plot(t_open,squeeze(Q_open(i_state,pix_num,:,:)),':c','Linewidth',2);
    h2 = plot(t_true,squeeze(Q_true(i_state,pix_num,:)),'-k','Linewidth',2);
    h3 = plot(t_open,squeeze(mean(Q_open(i_state,pix_num,:,:),4)),'-m','LineWidth',2);
    xlabel('DOY'); 
    legend([h2; h3; h1(1)],'Truth','OL Mean','OL Replicates')
    grid on
    if i_state==1
        ylabel('R_n [W m^{-2}]');
    elseif i_state==2
        ylabel('H [W m^{-2}]'); 
    else
        ylabel('LE [W m^{-2}]'); 
    end
    
    subplot(1,2,2)
    hold on
    h1 = plot(t,squeeze(Q(i_state,pix_num,:,:)),':c','LineWidth',2);
    h2 = plot(t_true,squeeze(Q_true(i_state,pix_num,:)),'-k','Linewidth',2);
    h3 = plot(t,squeeze(mean(Q(i_state,pix_num,:,:),4)),'-m','LineWidth',2);
    xlabel('DOY');
    legend([h2; h3; h1(1)],'Truth','EnKF Mean','EnKF Replicates')
    grid on

end

% Show some mapping results
if M>1
    for i_t=1:length(t_meas)
        ttrue_meas_index(i_t)=find(abs(t_true-t_meas(i_t))<(t_true(2)-t_true(1))/2);
        tmeas_index(i_t,:)=find(abs(t-t_meas(i_t))<(t(2)-t(1))/2);
    end

    % Plot soil moisture maps at measurement times only (for brevity)
    W1_lim=[0.25 0.6]; W1_range_lim=[0.0 0.5];
    W2_lim=[0.3 0.5]; W2_range_lim=[0.1 0.3];

    % Surface soil moisture maps
    for i=1:length(t_meas)
        figure(8+i);clf    
        subplot(2,4,1); imagesc(reshape(squeeze(Y_true(3,:,ttrue_meas_index(i))),N_lat,N_lon),W1_lim); axis image; colorbar; title(['True W_s at Meas. #: ' num2str(i)])
        subplot(2,4,5); imagesc(z_meas(i),W1_lim); axis image; colorbar; title('Measurement')
        subplot(2,4,2); imagesc(reshape(squeeze(mean(Y_open(3,:,ttrue_meas_index(i),:),4)),N_lat,N_lon),W1_lim); axis image; colorbar; title('OL Mean')
        subplot(2,4,6); imagesc(reshape(squeeze(max(Y_open(3,:,ttrue_meas_index(i),:),[],4) - min(Y_open(3,:,ttrue_meas_index(i),:),[],4)   ),N_lat,N_lon),W1_range_lim); axis image; colorbar; title('OL Range')
        subplot(2,4,3); imagesc(reshape(squeeze(mean(Y(3,:,tmeas_index(i,1),:),4)),N_lat,N_lon),W1_lim); axis image; colorbar; title('EnKF Prior Mean')
        subplot(2,4,7); imagesc(reshape(squeeze(max(Y(3,:,tmeas_index(i,1),:),[],4) - min(Y(3,:,tmeas_index(i,1),:),[],4)),N_lat,N_lon),W1_range_lim); axis image; colorbar; title('EnKF Prior Range')
        subplot(2,4,4); imagesc(reshape(squeeze(mean(Y(3,:,tmeas_index(i,2),:),4)),N_lat,N_lon),W1_lim); axis image; colorbar; title('EnKF Posterior Mean')
        subplot(2,4,8); imagesc(reshape(squeeze(max(Y(3,:,tmeas_index(i,2),:),[],4) - min(Y(3,:,tmeas_index(i,2),:),[],4)),N_lat,N_lon),W1_range_lim); axis image; colorbar; title('EnKF Posterior Range')

    end

    % Rootzone soil moisture maps
    for i=1:length(t_meas)
        figure(8+length(t_meas)+i);clf
        subplot(2,4,1); imagesc(reshape(squeeze(Y_true(4,:,ttrue_meas_index(i))),N_lat,N_lon),W2_lim); axis image; colorbar; title(['True W_d at Meas. #: ' num2str(i)])
        subplot(2,4,2); imagesc(reshape(squeeze(mean(Y_open(4,:,ttrue_meas_index(i),:),4)),N_lat,N_lon),W2_lim); axis image; colorbar; title('OL Mean')
        subplot(2,4,6); imagesc(reshape(squeeze(max(Y_open(4,:,ttrue_meas_index(i),:),[],4) - min(Y_open(4,:,ttrue_meas_index(i),:),[],4)   ),N_lat,N_lon),W2_range_lim); axis image; colorbar; title('OL Range')
        subplot(2,4,3); imagesc(reshape(squeeze(mean(Y(4,:,tmeas_index(i,1),:),4)),N_lat,N_lon),W2_lim); axis image; colorbar; title('EnKF Prior Mean')
        subplot(2,4,7); imagesc(reshape(squeeze(max(Y(4,:,tmeas_index(i,1),:),[],4) - min(Y(4,:,tmeas_index(i,1),:),[],4)),N_lat,N_lon),W2_range_lim); axis image; colorbar; title('EnKF Prior Range')
        subplot(2,4,4); imagesc(reshape(squeeze(mean(Y(4,:,tmeas_index(i,2),:),4)),N_lat,N_lon),W2_lim); axis image; colorbar; title('EnKF Posterior Mean')
        subplot(2,4,8); imagesc(reshape(squeeze(max(Y(4,:,tmeas_index(i,2),:),[],4) - min(Y(4,:,tmeas_index(i,2),:),[],4)),N_lat,N_lon),W2_range_lim); axis image; colorbar; title('EnKF Posterior Range')

    end
end

toc;
