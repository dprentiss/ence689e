%PS6
Prob2_single_state_truth_gen;
t = 0:150;
tz = z + 1
hold on;
plot(t, ytrue);
scatter(tmeas,z);

figure
plot(t, utrue)

noise = zeros(1,150);
for i = 1:150
    noise(i) = noise(i) + normrnd(0,w_variance);
end
noise = [0 noise];

j = 0;
for rho = -1:0.01:1
    j = j + 1;
    u = utrue*rho + noise;
    u = [0 u];
    rhos(j,:) = [rho, fct_bias(utrue, u(1:end-1)), fct_RMSE(utrue, u(1:end-1))]
end

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