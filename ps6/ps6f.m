Prob2_dual_state_truth_gen
alpha(4)=0.99
Prob2_open_loop

figure
hold on
plot(ttrue,ytrue(1,:))
plot(ttrue,squeeze(ybar_FULL(1,:)))
plot(ttrue, 2*sqrt(squeeze(Cyy_FULL(1,1,:))))
plot(ttrue, -2*sqrt(squeeze(Cyy_FULL(1,1,:))))

figure
hold on
plot(ttrue,utrue(1,:))
plot(ttrue,squeeze(ybar_FULL(2,:)))
plot(ttrue, 2*sqrt(squeeze(Cyy_FULL(2,2,:))))
plot(ttrue, -2*sqrt(squeeze(Cyy_FULL(2,2,:))))

fct_bias(squeeze(ybar_FULL(1,:)),ytrue(1,:))
fct_RMSE(squeeze(ybar_FULL(1,:)),ytrue(1,:))
fct_bias(squeeze(ybar_FULL(2,:)),utrue)
fct_RMSE(squeeze(ybar_FULL(2,:)),utrue)