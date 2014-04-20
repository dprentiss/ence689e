ubar = [0 0 70 70 0 0 0 0 0 0];
T = 0:9;
plot(T,ubar)
axis([0 9 -10 80])
title('Nominal precipitation forcing')
ylabel('mm water')
xlabel('time step')

rng(0);
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