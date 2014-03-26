function bias_return_var = fct_bias(observed,simulated)

% Function used to compute the bias between a series of simulated
% values and observed values computed as SIMULATED-OBSERVED.
%
% Inputs:
%   Observation vector: observed
%   Model output vector: simulated
%
% Author:   Bart Forman

% Capture and coerce into column vectors
observed = observed(:);
simulated = simulated(:);

% Determine coincident indices
ind_OK = find(~isnan(observed) & ~isnan(simulated));
if isempty(ind_OK)
    bias_return_var = NaN;
else
    % Compute residual using "valid" values only
    res = simulated(ind_OK) - observed(ind_OK);

    % Compute bias (set as the return value)
    bias_return_var = sum(res) / length(ind_OK);
end