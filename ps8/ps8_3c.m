%after running the true and EnKF scripts...

load true_output.mat imeas Y_true Q_true;
load EnKF_output.mat Y Q;
Ysize = fct_resize_y(imeas, mean(squeeze(Y),3));
Qsize = fct_resize_y(imeas, mean(squeeze(Q),3));
tmat = zeros(2,7);
for i = 1:4
    tmat(1,i) = fct_RMSE(squeeze(Y_true(i,:,:)), Ysize(i,:));
    tmat(2,i) = fct_bias(squeeze(Y_true(i,:,:)), Ysize(i,:));
end

for i = 5:7
    tmat(1,i) = fct_RMSE(squeeze(Q_true(i-4,:,:)), Qsize(i-4,:));
    tmat(2,i) = fct_bias(squeeze(Q_true(i-4,:,:)), Qsize(i-4,:));
end

digits(4);
latex(vpa(sym(tmat)))