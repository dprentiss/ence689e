function Y_est = fct_resize_y(imeas,Y)

% Function removes the prior estimates at the update times from
% the conditioned ensemble to facilitate the comparison against
% ensemble open-loop results.
%
%   Inputs: Temporal indices of measurement(s)
%           State vector (including prior estimates at measurement times)
%   Output: State vector (excluding prior estimates at measurement times)
%
% Author:   Bart Forman
% Created:  10 Feb 2012

Y_est = Y; % initialize the "full" matrix
for n = 1:length(imeas)
    Y_est(:,imeas(n),:) = []; % remove prior update times
end
