A = [0.8 0.1; 0.04 0.8];
G = [1; 0];
u = [0 0 70 70 0 0 0 0 0 0];
T = 0:10;

y_true = zeros(2,11);
y_true(:,1) = [50; 75];
y_b = zeros(2,11);
y_b(:,1)= [70; 80];

for i = 1:10
    y_true(:,i+1) = A*y_true(:,i) + G*u(i);
    y_b(:,i+1) = A*y_b(:,i) + G*u(i);
end

plot(T, y_true(1,:), '-sk', T, y_true(2,:), '-ok', ...
    T, y_b(1,:), ':sk', T, y_b(2,:), ':ok')

y1biasb = fct_bias(y_b(1,:),y_true(1,:))
y2biasb = fct_bias(y_b(2,:),y_true(2,:))
y1RMSEb = fct_RMSE(y_b(1,:),y_true(1,:))
y2RMSEb = fct_RMSE(y_b(2,:),y_true(2,:))

Cyy_c(:,:,1) = [400, 0; 0, 400];

for i = 1:10
    Cyy_c(:,:,i+1) = A*Cyy_c(:,:,i)*A';
end

y1biasc = fct_bias(y_b(1,:),y_true(1,:))
y2biasc = fct_bias(y_b(2,:),y_true(2,:))
y1RMSEc = fct_RMSE(y_b(1,:),y_true(1,:))
y2RMSEc = fct_RMSE(y_b(2,:),y_true(2,:))

yb1_up = y_b(1,:)' + 2*sqrt(squeeze(Cyy_c(1,1,:)));
yb1_lw = y_b(1,:)' - 2*sqrt(squeeze(Cyy_c(1,1,:)));
yb2_up = y_b(2,:)' + 2*sqrt(squeeze(Cyy_c(2,2,:)));
yb2_lw = y_b(2,:)' - 2*sqrt(squeeze(Cyy_c(2,2,:)));

plot(T, y_true(1,:), '-sk', T, y_true(2,:), '-ok', ...
    T, y_b(1,:), '-.sk', T, y_b(2,:), '-.sk', ...
    T, yb1_up, ':k', T, yb1_lw, ':k', ...
    T, yb2_up, ':k', T, yb2_lw, ':k')

rng(1);
y_E = zeros(2,1000,11);
y_E(:,:,1) = mvnrnd(y_b(:,1),Cyy_c(:,:,1),1000)';

for i = 1:10
    y_E(:,:,i+1) = A*y_E(:,:,i) + repmat(G*u(i),1,1000);
end

y_E_mean = mean(y_E(:,:,1),2)
y_E_covar = cov(y_E(:,:,1)')

% plot(T, y_true(1,:), '-k', T, y_true(2,:), '-k')
% hold on;
% for i = 1:100:1000
%      plot(T, squeeze(y_E(1,i,:)), ':k', T, squeeze(y_E(2,i,:)), ':k')
% end

