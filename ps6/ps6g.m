Prob2_dual_state_truth_gen
Prob2_KF

for kk = 1:length(imeas)
    ybar(:,imeas(kk)) = [];
end

fct_bias(squeeze(ybar(1,:)),ytrue(1,:))
fct_RMSE(squeeze(ybar(1,:)),ytrue(1,:))
fct_bias(squeeze(ybar(2,:)),utrue)
fct_RMSE(squeeze(ybar(2,:)),utrue)