
SWE_sum = 0;
for i = [2:744]
    SWE_sum = SWE_sum + SWE_maps(:,:,i) - SWE_maps(:,:,i-1);
end