% calculate r.v. x
m = [1; 2]; v = [0.09; 0.25];    % mean and var of x
mu = log((m.^2)./sqrt(v+m.^2));  % parameter mu
sigma = log(v./(m.^2)+1);        % parameter sigma
sigma = diag(sigma);             % create matrix from var
rng(1);                          % set random seed
w = mvnrnd(mu, sigma, 1000);     % calculate normal r.v. ensemble
x = exp(w);                      % calculate lognormal r.v. ensemble
mean(x), cov(x)                  % display mean and cov

% define F(x)
F = [1+0.25*x(:,1)-0.1*x(:,2).^2,...
    0.05*x(:,1)+0.1*x(:,1).^2+0.1*x(:,2).^2];
mean(F), cov(F)                  % display mean and cov
corrcoef(F)                      % display correlation