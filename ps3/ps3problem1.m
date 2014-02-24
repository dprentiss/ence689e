syms x;
pdf = sym(exp(-x));    % define exponential density function
cdf = int(pdf, 0, x);  % integrate pdf to find cdf
F = finverse(cdf);     % calculate inverse of cdf
rng(1);                % set random seed for reproducible results

plotNum = 0;
% calculate and plot ensembles generated with rand
for n = [10 100 1000 10000]
    bins = 0:0.1:10;
    r = rand(n, 1);                 % generate random probabilties
    mc = double(subs(F, 'x', r));   % calculate z values
    [count, loc] = hist(mc, bins);
    
    % plot figure 1
    figure(1)
    plotNum = plotNum + 1;
    subplot(2, 2, plotNum)
    bar(loc, count/n*10, 'EdgeColor', 'none', 'BarWidth', 1); hold on;
    d = 0:0.01:6;
    plot(d, exp(-d));
    axis([0, 6, 0, 1]);
    title(['exp(-x), n=', num2str(n)]);
    set(findall(gcf,'type','text'),'fontSize',18);
    set(gca,'fontSize',18);
end
mean(mc), var(mc)


plotNum = 0;
% calculate and plot ensembles generated with exprnd
for n = [10 100 1000 10000]
    bins = 0:0.1:10;
    mc = exprnd(1, n, 1);            % generate random probabilties
    [count, loc] = hist(mc, bins);
    
    % plot figure 2
    figure(2)
    plotNum = plotNum + 1;
    subplot(2, 2, plotNum)
    bar(loc, count/n*10, 'EdgeColor', 'none', 'BarWidth', 1); hold on;
    d = 0:0.01:6;
    plot(d, exp(-d));
    axis([0, 6, 0, 1]);
    title(['exp(-x), n=', num2str(n)]);
    set(findall(gcf,'type','text'),'fontSize',18);
    set(gca,'fontSize',18);
end
mean(mc), var(mc)