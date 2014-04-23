function RMSE_return_var = fct_RMSE(observed,simulated)

% Function used to compute the root mean squared error (RMSE) between
% a series of simulated values and observed values.
% 
% Inputs:
%   Observation vector: observed
%   Model output vector: simulated
%
% Author:   Bart Forman

% Capture and coerce into column vectors
observed = observed(:);
simulated = simulated(:);

% Determine coincident indices>0
ind_OK = find(~isnan(observed) & ~isnan(simulated));
if isempty(ind_OK)
    RMSE_return_var = NaN;
else
    % Compute residual using "valid" values only
    res = simulated(ind_OK) - observed(ind_OK);
    
    % Compute RMSE and set as the return value
    sum_res_squared = sum(res.*res);
    RMSE_return_var = sqrt(sum_res_squared/length(ind_OK));
end