SWE_sum = 0;
for i = [2:744]
    SWE_sum = SWE_sum + SWE_maps(:,:,i) - SWE_maps(:,:,i-1);
end
% generate labels from lat_deg and lon_deg
xTL = lon_deg(1:10:151);
xT = linspace(1, 151, numel(xTL));
set(gca, 'XTick', xT, 'XTickLabel', xTL)