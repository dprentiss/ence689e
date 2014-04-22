%nominal forcing
ubar = [0 0 70 70 0 0 0 0 0 0];
T = 0:9;
plot(T,ubar)
axis([0 9 -10 80])
title('Nominal precipitation forcing')
ylabel('mm water')
xlabel('time step')

rng(0);

%additive, uncorrelated error
u = zeros(10,10);
for i = 1:10
    for j = 1:10
        u(i,j) = ubar(j) + normrnd(0, 1);
    end
end

figure
plot(T,ubar)
hold on
for i = 1:10
    plot(T,u(i,:))
end
axis([0 9 -10 80])
title('Ensemble precipitation forcing')
ylabel('mm water')
xlabel('time step')

%additive, correlated error
u = zeros(10,10);
for i = 1:10
    u(i,:) = ubar(:) + normrnd(0, 1);
end

figure
plot(T,ubar)
hold on
for i = 1:10
    plot(T,u(i,:))
end
axis([0 9 -10 80])
title('Ensemble precipitation forcing')
ylabel('mm water')
xlabel('time step')

%multiplicative, uncorrelated error
for i = 1:10
    for j = 1:10
        u(i,j) = ubar(j) * lognrnd(1, 0.1);
    end
end

figure
plot(T,ubar)
hold on
for i = 1:10
    plot(T,u(i,:))
end
axis([0 9 -10 200])
title('Ensemble precipitation forcing')
ylabel('mm water')
xlabel('time step')

%additive, correlated error
u = zeros(10,10);
for i = 1:10
    u(i,:) = ubar(:) * lognrnd(0, 0.1);
end

figure
plot(T,ubar)
hold on
for i = 1:10
    plot(T,u(i,:))
end
axis([0 9 -10 90])
title('Ensemble precipitation forcing')
ylabel('mm water')
xlabel('time step')