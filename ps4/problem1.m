load('PS4_prob1');
A = x_meas.^2;
Y = z_meas;
alpha = inv(A'*A)*A'*Y

x_1 = linspace(min(x_meas(:,1)), max(x_meas(:,1)));
x_2 = linspace(min(x_meas(:,2)), max(x_meas(:,2)));
[X_1, X_2] = meshgrid(x_1,x_2);
Z=alpha(1).*X_1.^2 + alpha(2).*X_2.^2;
contourf(X_1, X_2, Z, 30);