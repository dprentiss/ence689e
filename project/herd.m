rng(0);
load('shoatsRFE.mat')

% create seasons of effective rainfall by random walk
numSeasons = 20;
RFEs = zeros(1,numSeasons);
RFEs(1) = 1;
m = .8;
v = 0.1;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
for i = 2:numSeasons
    RFEs(i) = RFEs(i-1) * lognrnd(mu, sigma);
end

% generate truth
herdEnsembleSize = 1;
h = zeros(numSeasons,4,herdEnsembleSize);
h(1,:,:,:,:) = [50; 50; 50; 50];
for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	h(i,1,:) = h(i-1,1,:) + h(i-1,3,:)...
   	 - round(interp1(RFE,salesFemale,RFEs(i-1))*h(i-1,1,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,1,:));

    % NBt+1 = conception rate(AFt)
    h(i,2,:) = round(interp1(RFE,conceptions,RFEs(i-1))*h(i-1,1,:));
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    h(i,3,:) = round(0.5*(h(i-1,2,:)...
     - round(interp1(RFE,mortImm,RFEs(i-1))*h(i-1,2,:))))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,3,:));
 

    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    h(i,4,:) = h(i-1,4,:) + h(i,3,:)...
     - round(interp1(RFE,salesMale,RFEs(i-1))*h(i-1,4,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,4,:));
end
herdSize = sum(h,2);
figure;
hold on;
for i = 1:herdEnsembleSize
    plot(herdSize(:,1,i));
end
figure;
plot(RFEs);
figure;
hold on;
clear hratio;
hratio(:,1,:) = h(:,1,:)./herdSize;
hratio(:,2,:) = h(:,2,:)./herdSize;
hratio(:,3,:) = h(:,3,:)./herdSize;
hratio(:,4,:) = h(:,4,:)./herdSize;
plot(squeeze(hratio(:,1,:)))
plot(squeeze(hratio(:,2,:)))
plot(squeeze(hratio(:,3,:)))
plot(squeeze(hratio(:,4,:)))




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
hAF = randi([herdMinAF herdMaxAF], 1, herdEnsembleSize);
hNB = randi([herdMinNB herdMaxNB], 1, herdEnsembleSize);
hYF = randi([herdMinYF herdMaxYF], 1, herdEnsembleSize);
hYM = randi([herdMinYM herdMaxYM], 1, herdEnsembleSize);

h = zeros(numSeasons,4,herdEnsembleSize);
h(1,:,:,:,:)=[hAF; hNB; hYF; hYM];

for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	h(i,1,:) = h(i-1,1,:) + h(i-1,3,:)...
   	 - round(interp1(RFE,salesFemale,RFEs(i-1))*h(i-1,1,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,1,:));

    % NBt+1 = conception rate(AFt)
    h(i,2,:) = round(interp1(RFE,conceptions,RFEs(i-1))*h(i-1,1,:));
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    h(i,3,:) = round(0.5*(h(i-1,2,:)...
     - round(interp1(RFE,mortImm,RFEs(i-1))*h(i-1,2,:))))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,3,:));
 

    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    h(i,4,:) = h(i-1,4,:) + h(i,3,:)...
     - round(interp1(RFE,salesMale,RFEs(i-1))*h(i-1,4,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,4,:));
end
herdSize = sum(h,2);
figure;
hold on;
for i = 1:herdEnsembleSize
    plot(herdSize(:,1,i));
end
figure;
plot(RFEs);
figure;
hold on;
clear hratio;
hratio(:,1,:) = h(:,1,:)./herdSize;
hratio(:,2,:) = h(:,2,:)./herdSize;
hratio(:,3,:) = h(:,3,:)./herdSize;
hratio(:,4,:) = h(:,4,:)./herdSize;
plot(squeeze(hratio(:,1,:)))
plot(squeeze(hratio(:,2,:)))
plot(squeeze(hratio(:,3,:)))
plot(squeeze(hratio(:,4,:)))