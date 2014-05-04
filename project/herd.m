clear all
rng(0);
load('shoatsRFE.mat')

% create seasons of effective rainfall forcing
numSeasons = 40;
RFEs = zeros(1,numSeasons);
RFEs(1) = 1;
m = .7;
v = 0.1;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
for i = 2:numSeasons
    RFEs(i) = RFEs(i-1) * lognrnd(mu, sigma);
end

% generate truth
herdEnsembleSize = 1;
y = zeros(numSeasons,4,herdEnsembleSize);
y(1,:,:,:,:) = [100; 100; 100; 100];
for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	y(i,1,:) = y(i-1,1,:) + y(i-1,3,:)...
   	 - round(interp1(RFE,salesFemale,RFEs(i-1))*y(i-1,1,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*y(i-1,1,:));

    % NBt+1 = conception rate(AFt)
    y(i,2,:) = round(interp1(RFE,conceptions,RFEs(i-1))*y(i-1,1,:));
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    y(i,3,:) = round(0.5*(y(i-1,2,:)...
     - round(interp1(RFE,mortImm,RFEs(i-1))*y(i-1,2,:))))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*y(i-1,3,:));
 
    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    y(i,4,:) = y(i-1,4,:) + y(i,3,:)...
     - round(interp1(RFE,salesMale,RFEs(i-1))*y(i-1,4,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*y(i-1,4,:));
end

ytrue = y;

% plot true total herd size
herdSize = sum(y,2);
figure(1);
hold on;
for i = 1:herdEnsembleSize
    plot(herdSize(:,1,i));
end
title('True Herd Size')

%plot forcing
figure(2)
plot(RFEs)
title('Effective Rainfall')

% plot true demogaphic ratios
% figure(3)
% clf(3)
% hold on
% clear yratio
% yratio(:,1,:) = y(:,1,:)./herdSize;
% yratio(:,2,:) = y(:,2,:)./herdSize;
% yratio(:,3,:) = y(:,3,:)./herdSize;
% yratio(:,4,:) = y(:,4,:)./herdSize;
% plot(squeeze(yratio(:,1,:)))
% plot(squeeze(yratio(:,2,:)))
% plot(squeeze(yratio(:,3,:)))
% plot(squeeze(yratio(:,4,:)))

% generate synthetic measurements of newborns
m = 1;
v = 0.01;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
measStep = 4;
numMeas = (numSeasons - mod(numSeasons, measStep))/measStep;
tmeas = zeros(1, numMeas);
zNB = zeros(1, numMeas);
for i = 1:numMeas
    tmeas(i) = i*measStep;
    zNB(i) = round(ytrue(tmeas(i),2,1) * lognrnd(mu, sigma));
end

% plot truth with measurements
figure(4)
clf(4)
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

% initialize herd ensemble with uniformly distibuted,
% uncorrelated herd demographic groups
% AF -- adult female
% NB -- newborn
% YF -- young female
% YM -- young male
% NOTE: Adult males, if present, are for studding only and are not modeled
herdEnsembleSize = 100;
herdMinAF = 0;
herdMaxAF = 100;
herdMinNB = 0;
herdMaxNB = 100;
herdMinYF = 0;
herdMaxYF = 100;
herdMinYM = 0;
herdMaxYM = 100;
yAF = randi([herdMinAF herdMaxAF], 1, herdEnsembleSize);
yNB = randi([herdMinNB herdMaxNB], 1, herdEnsembleSize);
yYF = randi([herdMinYF herdMaxYF], 1, herdEnsembleSize);
yYM = randi([herdMinYM herdMaxYM], 1, herdEnsembleSize);

y = zeros(numSeasons,4,herdEnsembleSize);
y(1,:,:,:,:)=[yAF; yNB; yYF; yYM];

m = 1;
v = 0.001;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	y(i,1,:) = y(i-1,1,:) + y(i-1,3,:)...
   	 - round(interp1(RFE,salesFemale,RFEs(i-1))*y(i-1,1,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*y(i-1,1,:));
 
    y(i,1,:) = round(squeeze(y(i,1,:)) .* lognrnd(mu, sigma, herdEnsembleSize, 1));

    % NBt+1 = conception rate(AFt)
    y(i,2,:) = round(interp1(RFE,conceptions,RFEs(i-1))*y(i-1,1,:));
    
    y(i,2,:) = round(squeeze(y(i,2,:)) .* lognrnd(mu, sigma, herdEnsembleSize, 1));
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    y(i,3,:) = round(0.5*(y(i-1,2,:)...
     - round(interp1(RFE,mortImm,RFEs(i-1))*y(i-1,2,:))))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*y(i-1,3,:));
    
    y(i,3,:) = round(squeeze(y(i,3,:)) .* lognrnd(mu, sigma, herdEnsembleSize, 1));

    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    y(i,4,:) = y(i-1,4,:) + y(i,3,:)...
     - round(interp1(RFE,salesMale,RFEs(i-1))*y(i-1,4,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*y(i-1,4,:));
 
    y(i,4,:) = round(squeeze(y(i,4,:)) .* lognrnd(mu, sigma, herdEnsembleSize, 1));
 
end

% plot openloop
herdSize = sum(y,2);
figure(5);
clf(5)
hold on;
for i = 1:herdEnsembleSize
    plot(herdSize(:,1,i));
end
% figure;
% hold on;
% clear yratio;
% yratio(:,1,:) = y(:,1,:)./herdSize;
% yratio(:,2,:) = y(:,2,:)./herdSize;
% yratio(:,3,:) = y(:,3,:)./herdSize;
% yratio(:,4,:) = y(:,4,:)./herdSize;
% plot(squeeze(yratio(:,1,:)))
% plot(squeeze(yratio(:,2,:)))
% plot(squeeze(yratio(:,3,:)))
% plot(squeeze(yratio(:,4,:)))

% Ensemble Kalman update

% Multiplicative error mean and covariance
m = 1;
v = 0.1;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

% K = Cyz[Czz+Cvv]^(-1)

% Cyz

%ybar = repmat(mean(y(i,:,:),3)',1,100);
%squeeze(y(i,:,:)) - ybar;
%zpNB = squeeze(y(i,2,:))
%zactNB = round(zNB(i) * lognrnd(mu, sigma, herdEnsembleSize, 1));