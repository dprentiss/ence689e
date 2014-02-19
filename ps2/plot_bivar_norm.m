% This function plots a bivariate normal pdf.
% 
% Inputs:
%   Mean vector: xbar
%   Covariance Matrix: Cxx
%
% Author:   Bart Forman

function plot_bivar_norm(xbar,Cxx)

% Compute mean values
mu1 = xbar(1);
mu2 = xbar(2);

% Compute standard deviations
sigma1 = sqrt(Cxx(1,1));
sigma2 = sqrt(Cxx(2,2));

% Compute the correlation coefficient
rho = Cxx(2,1) / sigma1 / sigma2;

% Generate a range of values for x1, x2
x1 = linspace(mu1-3*sigma1,mu1+3*sigma1,100);
x2 = linspace(mu2-3*sigma2,mu2+3*sigma2,100);

% Compute z at each point within the probability space of interest
for i=1:length(x1)
    z(i,:) = ... % NOTE: "..." simply means continue on the next line
        (x1(i)-mu1).^2/sigma1^2 - 2*rho*(x1(i)-mu1) .* ...
        (x2-mu2)/sigma1/sigma2 + (x2-mu2).^2/sigma2^2;
end

% Compute the bivariate probability distribution
p_unnormalized = 1/(2*pi*sigma1*sigma2*sqrt(1-rho^2)) * exp(-z/(2*(1-rho^2)));

% Normalize p so that the integrated pdf equals 1.0
p = p_unnormalized / sum(p_unnormalized(:));

% Generate desired plots
subplot(1,2,1)
mesh(x2,x1,p)
axis square % make x- and y-axis units equivalent
xlim([min(x2) max(x2)]) % set y-axis limits to the extent of x2
ylim([min(x1) max(x1)]) % set x-axis limits to the extent of x1
xlabel('X_2')
ylabel('X_1')
zlabel('Probability [-]')

subplot(1,2,2)
contour(x2,x1,p)
axis square % make x- and y-axis units equivalent
xlim([min(x2) max(x2)]) % set y-axis limits to the extent of x2
ylim([min(x1) max(x1)]) % set x-axis limits to the extent of x1
xlabel('X_2')
ylabel('X_1')
