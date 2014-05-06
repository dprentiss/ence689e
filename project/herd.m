% AF -- adult female
% NB -- newborn
% YF -- young female
% YM -- young male
% NOTE: Adult males, if present, are for studding only and are not modeled
clear all
rng(0);
load('shoatsRFE.mat')
figNum = 0;

mortRate = 0.001;

% Model error
m = 1;
v = 0.001;
ymu = log((m^2)/sqrt(v+m^2));
ysigma = sqrt(log(v/(m^2)+1));

% measurement error
m = 1;
v = 0.01;
Cvv = v;
zmu = log((m^2)/sqrt(v+m^2));
zsigma = sqrt(log(v/(m^2)+1));

% effective rainfall forcing, synthetic truth
numSeasons = 40;
RFEs = zeros(1,numSeasons);
RFEs(1) = 1;
m = .8;
%m = 0.2
v = 0.1;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
for i = 2:numSeasons
    %RFEs(i) = 1;
    RFEs(i) = RFEs(i-1) * lognrnd(mu, sigma);
    %RFEs(i) = lognrnd(mu, sigma);
end

% generate truth
herdEnsembleSize = 1;
y = zeros(numSeasons,4,herdEnsembleSize);
y(1,:,:,:,:) = [100; 100; 100; 100];
for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	y(i,1,:) = y(i-1,1,:) + y(i-1,3,:)...
   	 - interp1(RFE,salesFemale,RFEs(i-1))*y(i-1,1,:)...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,1,:);
     %- (sum(y(i-1,:,:),2)/1000).*interp1(RFE,mortMat,RFEs(i-1)).*y(i-1,1,:);

    % NBt+1 = conception rate(AFt)
    y(i,2,:) = interp1(RFE,conceptions,RFEs(i-1))*y(i-1,1,:);
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    y(i,3,:) = 0.5*(y(i-1,2,:)...
     - interp1(RFE,mortImm,RFEs(i-1))*y(i-1,2,:))...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,3,:);
 
    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    y(i,4,:) = y(i-1,4,:) + y(i,3,:)...
     - interp1(RFE,salesMale,RFEs(i-1))*y(i-1,4,:)...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,4,:);
end
ytrue = y;

% plot true total herd size
herdSize = sum(y,2);
% figNum += 1;
% figure(figNum);
% clg(figNum);
% hold on;
% for i = 1:herdEnsembleSize
%     plot(herdSize(:,1,i));
% end
% title('True Herd Size')

%plot forcing
figNum = figNum + 1;
figure(figNum)
clf(figNum)
plot(RFEs)
title('Effective Rainfall')

% plot true demogaphic ratios
figNum = figNum + 1;
figure(figNum)
clf(figNum)
hold on
clear yratio
yratio(:,1,:) = y(:,1,:)./herdSize;
yratio(:,2,:) = y(:,2,:)./herdSize;
yratio(:,3,:) = y(:,3,:)./herdSize;
yratio(:,4,:) = y(:,4,:)./herdSize;
plot(squeeze(yratio(:,1,:)))
plot(squeeze(yratio(:,2,:)))
plot(squeeze(yratio(:,3,:)))
plot(squeeze(yratio(:,4,:)))

% generate synthetic measurements of newborns
measStep = 4;
numMeas = (numSeasons - mod(numSeasons, measStep))/measStep;
tmeas = zeros(1, numMeas);
zNB = zeros(1, numMeas);
for i = 1:numMeas
    tmeas(i) = i*measStep;
    zNB(i) = ytrue(tmeas(i),2,1) * lognrnd(zmu, zsigma);
end

% plot truth with measurements
figNum = figNum + 1;
figure(figNum)
clf(figNum)
subplot(4,1,1), plot(ytrue(:,1))
title('Adult Females')
subplot(4,1,2), plot(ytrue(:,2))
title('Newborns and Measurements')
hold on
subplot(4,1,2), plot(tmeas, zNB, 'o')
subplot(4,1,3), plot(ytrue(:,3))
title('Young Females')
subplot(4,1,4), plot(ytrue(:,4))
title('Young Males')

% % general dynamics
% RFEg = ones(1,numSeasons);
% for i = 0:0.1:1
%     RFEg = ones(1,numSeasons);
%     RFEg = RFEg * i;
%     herdEnsembleSize = 100;
%     herdMinAF = 0;
%     herdMaxAF = 200;
%     herdMinNB = 0;
%     herdMaxNB = 200;
%     herdMinYF = 0;
%     herdMaxYF = 200;
%     herdMinYM = 0;
%     herdMaxYM = 200;
%     yAF = randi([herdMinAF herdMaxAF], 1, herdEnsembleSize);
%     yNB = randi([herdMinNB herdMaxNB], 1, herdEnsembleSize);
%     yYF = randi([herdMinYF herdMaxYF], 1, herdEnsembleSize);
%     yYM = randi([herdMinYM herdMaxYM], 1, herdEnsembleSize);
%     
%     y = zeros(numSeasons,4,herdEnsembleSize);
%     y(1,:,:,:,:)=[yAF; yNB; yYF; yYM];
%     
%     
%     for k = 2:numSeasons % for every seasons
%         
%         % AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
%         y(k,1,:) = y(k-1,1,:) + y(k-1,3,:)...
%             - interp1(RFE,salesFemale,RFEg(k-1))*y(k-1,1,:)...
%             - interp1(RFE,mortMat,RFEg(k-1))*y(k-1,1,:);
%         
%         y(k,1,:) = squeeze(y(k,1,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
%         
%         % NBt+1 = conception rate(AFt)
%         y(k,2,:) = interp1(RFE,conceptions,RFEg(k-1))*y(k-1,1,:);
%         
%         y(k,2,:) = squeeze(y(k,2,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
%         
%         % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
%         y(k,3,:) = 0.5*(y(k-1,2,:)...
%             - interp1(RFE,mortImm,RFEg(k-1))*y(k-1,2,:))...
%             - interp1(RFE,mortMat,RFEg(k-1))*y(k-1,3,:);
%         
%         y(k,3,:) = squeeze(y(k,3,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
%         
%         % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
%         y(k,4,:) = y(k-1,4,:) + y(k,3,:)...
%             - interp1(RFE,salesMale,RFEg(k-1))*y(k-1,4,:)...
%             - interp1(RFE,mortMat,RFEg(k-1))*y(k-1,4,:);
%         
%         y(k,4,:) = squeeze(y(k,4,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
%         
%     end
%        
%     % plot openloop
%     herdSize = sum(y,2);
%     figNum = figNum + 1;
%     figure(figNum)
%     clf(figNum)
%     subplot(4,1,1)
%     hold on
%     for j = 1:herdEnsembleSize
%         plot(y(:,1,j),'c:');
%     end
%     plot(ytrue(:,1))
%     title(sprintf('Adult Females, RFE = %f',i))
%     
%     subplot(4,1,2)
%     hold on
%     for j = 1:herdEnsembleSize
%         plot(y(:,2,j),'c:');
%     end
%     plot(ytrue(:,2))
%     title(sprintf('Newborns and Measurements, RFE = %f',i))
%     subplot(4,1,2), plot(tmeas, zNB, 'o')
%     
%     subplot(4,1,3)
%     title(sprintf('Young Females, RFE = %f',i))
%     hold on
%     for j = 1:herdEnsembleSize
%         plot(y(:,3,j),'c:');
%     end
%     plot(ytrue(:,3))
%     
%     subplot(4,1,4)
%     hold on
%     for j = 1:herdEnsembleSize
%         plot(y(:,4,j),'c:');
%     end
%     plot(ytrue(:,4))
%     title(sprintf('Young Males, RFE = %f',i))
%     
%     figNum = figNum + 1;
%     figure(figNum)
%     clf(figNum)
%     hold on
%     clear yratio
%     yratio(:,1,:) = y(:,1,:)./herdSize;
%     yratio(:,2,:) = y(:,2,:)./herdSize;
%     yratio(:,3,:) = y(:,3,:)./herdSize;
%     yratio(:,4,:) = y(:,4,:)./herdSize;
%     plot(squeeze(yratio(:,1,:)))
%     plot(squeeze(yratio(:,2,:)))
%     plot(squeeze(yratio(:,3,:)))
%     plot(squeeze(yratio(:,4,:)))
% end

% initialize herd ensemble with uniformly distibuted,
% uncorrelated herd demographic groups
herdEnsembleSize = 100;
herdMinAF = 100;
herdMaxAF = 200;
herdMinNB = 100;
herdMaxNB = 200;
herdMinYF = 100;
herdMaxYF = 200;
herdMinYM = 100;
herdMaxYM = 200;
yAF = randi([herdMinAF herdMaxAF], 1, herdEnsembleSize);
yNB = randi([herdMinNB herdMaxNB], 1, herdEnsembleSize);
yYF = randi([herdMinYF herdMaxYF], 1, herdEnsembleSize);
yYM = randi([herdMinYM herdMaxYM], 1, herdEnsembleSize);

y = zeros(numSeasons,4,herdEnsembleSize);
y(1,:,:,:,:)=[yAF; yNB; yYF; yYM];

for i = 2:numSeasons % for every seasons 
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	y(i,1,:) = y(i-1,1,:) + y(i-1,3,:)...
   	 - interp1(RFE,salesFemale,RFEs(i-1))*y(i-1,1,:)...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,1,:);
 
    y(i,1,:) = squeeze(y(i,1,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);

    % NBt+1 = conception rate(AFt)
    y(i,2,:) = interp1(RFE,conceptions,RFEs(i-1))*y(i-1,1,:);
    
    y(i,2,:) = squeeze(y(i,2,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    y(i,3,:) = 0.5*(y(i-1,2,:)...
     - interp1(RFE,mortImm,RFEs(i-1))*y(i-1,2,:))...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,3,:);
    
    y(i,3,:) = squeeze(y(i,3,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);

    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    y(i,4,:) = y(i-1,4,:) + y(i,3,:)...
     - interp1(RFE,salesMale,RFEs(i-1))*y(i-1,4,:)...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,4,:);
 
    y(i,4,:) = squeeze(y(i,4,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
 
end

yOLmean = mean(y,3);

% plot openloop
herdSize = sum(y,2);
figNum = figNum + 1;
figure(figNum)
clf(figNum)
subplot(4,1,1)
hold on
for i = 1:herdEnsembleSize
    plot(y(:,1,i),'c:')
end
plot(ytrue(:,1),'k-')
plot(yOLmean(:,1))
title('Adult Females')

subplot(4,1,2)
hold on
for i = 1:herdEnsembleSize
    plot(y(:,2,i),'c:');
end
plot(ytrue(:,2),'k-')
plot(yOLmean(:,2))
title('Newborns and Measurements')
subplot(4,1,2), plot(tmeas, zNB, 'o')

subplot(4,1,3)
title('Young Females')
hold on
for i = 1:herdEnsembleSize
    plot(y(:,3,i),'c:');
end
plot(yOLmean(:,3))
plot(ytrue(:,3),'k-')

subplot(4,1,4)
hold on
for i = 1:herdEnsembleSize
    plot(y(:,4,i),'c:');
end
plot(ytrue(:,4),'k-')
plot(yOLmean(:,4))
title('Young Males')

figNum = figNum + 1;
figure(figNum)
clf(figNum)
hold on
clear yratio
yratio(:,1,:) = y(:,1,:)./herdSize;
yratio(:,2,:) = y(:,2,:)./herdSize;
yratio(:,3,:) = y(:,3,:)./herdSize;
yratio(:,4,:) = y(:,4,:)./herdSize;
plot(squeeze(yratio(:,1,:)))
plot(squeeze(yratio(:,2,:)))
plot(squeeze(yratio(:,3,:)))
plot(squeeze(yratio(:,4,:)))

%%%%% Ensemble Kalman Filter

y = zeros(numSeasons,4,herdEnsembleSize);
y(1,:,:,:,:)=[yAF; yNB; yYF; yYM];

for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	y(i,1,:) = y(i-1,1,:) + y(i-1,3,:)...
   	 - interp1(RFE,salesFemale,RFEs(i-1))*y(i-1,1,:)...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,1,:);
 
    y(i,1,:) = squeeze(y(i,1,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);

    % NBt+1 = conception rate(AFt)
    y(i,2,:) = interp1(RFE,conceptions,RFEs(i-1))*y(i-1,1,:);
    
    y(i,2,:) = squeeze(y(i,2,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    y(i,3,:) = 0.5*(y(i-1,2,:)...
     - interp1(RFE,mortImm,RFEs(i-1))*y(i-1,2,:))...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,3,:);
    
    y(i,3,:) = squeeze(y(i,3,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);

    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    y(i,4,:) = y(i-1,4,:) + y(i,3,:)...
     - interp1(RFE,salesMale,RFEs(i-1))*y(i-1,4,:)...
     - interp1(RFE,mortMat,RFEs(i-1))*y(i-1,4,:);
 
    y(i,4,:) = squeeze(y(i,4,:)) .* lognrnd(ymu, ysigma, herdEnsembleSize, 1);
    
    % Kalman update
    if find(tmeas == i)
        % K = Cyz[Czz+Cvv]^(-1)
        ybar = repmat(mean(y(i,:,:),3)',1,herdEnsembleSize);
        zp = squeeze(y(i,2,:));
        zbar = repmat(mean(zp),herdEnsembleSize,1);
        Cyz = (squeeze(y(i,:,:)) - ybar) * (zp - zbar) / (herdEnsembleSize - 1);
        Czz = (zp - zbar)'*(zp - zbar) / (herdEnsembleSize - 1);
        K = Cyz/(Czz + Cvv);
        
        zact = zNB(i/measStep) * lognrnd(zmu, zsigma, herdEnsembleSize, 1);
        y(i,:,:) = y(i,:,:) + permute(K*(zact-zp)', [3 1 2]);
    end
end

yKFmean = mean(y,3);

% plot data assimilation
herdSize = sum(y,2);
figNum = figNum + 1;
figure(figNum)
clf(figNum)
subplot(4,1,1)
hold on
for i = 1:herdEnsembleSize
    plot(y(:,1,i),'c:')
end
plot(ytrue(:,1),'k-')
plot(yKFmean(:,1))
title('Adult Females')

subplot(4,1,2)
hold on
for i = 1:herdEnsembleSize
    plot(y(:,2,i),'c:');
end
plot(ytrue(:,2),'k-')
plot(yKFmean(:,2))
title('Newborns and Measurements')
subplot(4,1,2), plot(tmeas, zNB, 'o')

subplot(4,1,3)
title('Young Females')
hold on
for i = 1:herdEnsembleSize
    plot(y(:,3,i),'c:');
end
plot(yKFmean(:,3))
plot(ytrue(:,3),'k-')

subplot(4,1,4)
hold on
for i = 1:herdEnsembleSize
    plot(y(:,4,i),'c:');
end
plot(ytrue(:,4),'k-')
plot(yKFmean(:,4))
title('Young Males')

fct_bias(ytrue, yKFmean)
fct_bias(ytrue, yOLmean)
fct_RMSE(ytrue, yKFmean)
fct_RMSE(ytrue, yOLmean)

hold off



