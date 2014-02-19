%-- 2/18/2014 8:07 PM --%
yt = [10; 15];
Cyy = [9 10; 10 25];
A = [.9 -.1;.2 .75];
G = eye(2,2);
u = [3; 1];
Cuu = [2 0; 0 1];
b = [0.50; 0.75];
yt1 = A*yt + G*u + b;
Cyy1 = A*Cyy*A' + G*Cuu*G';
plot_bivar_norm(yt1,Cyy1);
sd1 = sqrt(Cyy1(1,1));        % standard deviation yt1
sd2 = sqrt(Cyy1(2,2));        % standard deviation yt2
rho1 = Cyy1(1,2)/(sd1 * sd2); % correlation coefficient