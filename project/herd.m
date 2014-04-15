rng(1);

% create seasons of effective rainfall by random walk
numSeasons
RFE = zeros(1,100);
RFE(1) = 0;
for i = 2:100
    RFE(i) = RFE(i-1) + normrnd(0,0.05);
end

% initialize herd ensemble with uniformly distibuted,
% uncorrelated herd demographic groups
% AF -- adult female
% NB -- newborn
% YF -- young female
% YM -- young female
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
h(1,:,:,:,:)=[hAF; hNB; hYF; hYM];

% propagate herd forward model
for i