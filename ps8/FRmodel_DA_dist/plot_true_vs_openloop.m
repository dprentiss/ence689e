% Script use to plot open-loop results relative to the "truth".
%
% Author:   Bart Forman
% Created:  10 Feb 2012

clear;  tic;

load true_output.mat
load prob_open_loop_output.mat

[N,M,P] = size(Y_true);

% Prompt user for input
pix_num = input(['Enter pixel number for time-series',...
    ' plotting (1...' num2str(M) '): ']);

%% Time-series plots
figure(1); clf

% Plot model state
for i=1:4
    subplot(2,4,i)
    hold on
    h1 = plot(t_open,squeeze(Y_open(i,pix_num,:,:)),':c','LineWidth',2);
    h2 = plot(t_true,squeeze(Y_true(i,pix_num,:)),'-k','LineWidth',2);
    h3 = plot(t_open,squeeze(mean(Y_open(i,pix_num,:,:),4)),'LineWidth',2);
    grid on
    xlabel('DOY','FontSize',12)
    if i==1
        ylabel('T_s [K]','FontSize',12); title(['Pixel #: ' num2str(pix_num)],'FontSize',12)
    elseif i==2
        ylabel('T_d [K]','FontSize',12); 
    elseif i==3
        ylabel('W_s [-]','FontSize',12); 
    else
        ylabel('W_d [-]','FontSize',12);
        legend([h2; h3; h1(1)],'Truth','OL Mean','OL Replicate')
    end
end

% Plot model fluxes
for i=1:3
    if i<4
        subplot(2,4,i+4)
        hold on
        h1 = plot(t_open,squeeze(Q_open(i,pix_num,:,:)),':c','LineWidth',2);
        h2 = plot(t_true,squeeze(Q_true(i,pix_num,:)),'-k','LineWidth',2);
        h3 = plot(t_open,squeeze(mean(Q_open(i,pix_num,:,:),4)),'LineWidth',2);
        grid on
    else
        G_true = Q_true(1,:,:)-Q_true(2,:,:)-Q_true(3,:,:);
        G_open = Q_open(1,:,:,:)-Q_open(2,:,:,:)-Q_open(3,:,:,:);
        subplot(2,4,i+4)
        hold on
        h1 = plot(t_open,squeeze(G_open(1,pix_num,:,:)),':c','LineWidth',2);
        h2 = plot(t_true,squeeze(G_true(1,pix_num,:)),'-k','LineWidth',2);
        h3 = plot(t_open,squeeze(mean(G_open(1,pix_num,:,:),4)),'-b','LineWidth',2);
        axis([sim_Day_beg sim_Day_end -50 50])
    end
    xlabel('DOY','FontSize',12)
    if i==1
        ylabel('R_n [W m^{-2}]','FontSize',12);
    elseif i==2
        ylabel('H [W m^{-2}]','FontSize',12); 
    elseif i==3
        ylabel('LE [W m^{-2}]','FontSize',12); 
    else
        ylabel('G [W m^{-2}]','FontSize',12);
    end
end

%% Spatial maps
if M>1  % Only plot maps if more than one pixel
    
    % Map plots
    figure(2);clf
    Yt_max=max(max(Y_true,[],3),[],2);
    Yt_min=min(min(Y_true,[],3),[],2);
    Yo_max=max(max(max(Y_open,[],4),[],3),[],2);
    Yo_min=min(min(min(Y_open,[],4),[],3),[],2);
    Y_max=max([Yt_max Yo_max],[],2);
    Y_min=min([Yt_min Yo_min],[],2);
    
    for tt=1:10:length(t_open)
    figure(2)
        for i=1:4
            subplot(2,4,i)
            imagesc(reshape(squeeze(Y_true(i,:,tt)),N_lat,N_lon),[Y_min(i) Y_max(i)]);
            colorbar; axis image
            if i==1
                title('T_s [K]','FontSize',12)
                ylabel(['True' '( t = ' num2str(floor(10*t_open(tt))/10) ')'],'FontSize',12)
            elseif i==2
                title('T_d [K]','FontSize',12)
            elseif i==3
                title('W_s [-]','FontSize',12)
            else
                title('W_d [-]','FontSize',12)
            end
        end

        for i=1:4
            subplot(2,4,i+4)
            imagesc(reshape(squeeze(mean(Y_open(i,:,tt,:),4)),N_lat,N_lon),[Y_min(i) Y_max(i)]);
            colorbar
            axis image
            if i==1
                ylabel(['Open-Loop Mean' '( t = ' num2str(floor(10*t_open(tt))/10) ')'],'FontSize',12)
            end
        end

        figure(1)
        for i=1:4
            subplot(2,4,i)
            plot(t_true(tt),squeeze(Y_true(i,pix_num,tt)),'+')
            hold on;
        end

        %pause(1)

    end
end

toc
