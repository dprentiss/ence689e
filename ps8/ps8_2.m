EnKF_demo
yest = fct_resize_y(imeas,Y);
yest = mean(yest, 3);
fct_RMSE([yest(1,:);yest(2,:)],ytrue)
fct_bias([yest(1,:);yest(2,:)],ytrue)

EnKF_demo_2
yest = fct_resize_y(imeas,Y);
yest = mean(yest, 3);
fct_RMSE([yest(1,:);yest(2,:)],ytrue)
fct_bias([yest(1,:);yest(2,:)],ytrue)

EnKF_demo_ol
yest = fct_resize_y(imeas,Y);
yest = mean(yest, 3);
fct_RMSE([yest(1,:);yest(2,:)],ytrue)
fct_bias([yest(1,:);yest(2,:)],ytrue)