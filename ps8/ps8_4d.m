%after running the true and EnKF scripts...


load true_output.mat imeas Y_true Q_true;
load EnKF_output.mat Y Q;
tmat = zeros(8,7);

for p = 1:4
    Ysize = fct_resize_y(imeas, mean(squeeze(Y(:,p,:,:)),3));
    Qsize = fct_resize_y(imeas, mean(squeeze(Q(:,p,:,:)),3));

    for i = 1:4
        tmat(2*p-1,i) = fct_RMSE(squeeze(Y_true(i,p,:)), Ysize(i,:));
        tmat(2*p,i) = fct_bias(squeeze(Y_true(i,p,:)), Ysize(i,:));
    end
    
    for i = 5:7
        tmat(2*p-1,i) = fct_RMSE(squeeze(Q_true(i-4,p,:)), Qsize(i-4,:));
        tmat(2*p,i) = fct_bias(squeeze(Q_true(i-4,p,:)), Qsize(i-4,:));
    end
end

digits(4);
latex(vpa(sym(tmat)))