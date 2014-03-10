A = [0.8 0.1; 0.04 0.8];
G = [1; 0];
u = [0 0 70 70 0 0 0 0 0 0];
T = 1:11;

y_a = cell(1,11);
y_a{1} = [50; 75];
y_b = cell(1,11);
y_b{1}= [70; 80];

for i = 1:10
    y_a{i+1} = A*y_a{i} + G*u(i);
    y_b{i+1} = A*y_b{i} + G*u(i);
end

y_a_mat = cell2mat(y_a);
y_true_1 = y_a_mat(1,:);
y_true_2 = y_a_mat(2,:);

y_b_mat = cell2mat(y_b);
y_b_1 = y_b_mat(1,:);
y_b_2 = y_b_mat(2,:);

plot(T, y_true_1, 'black', T, y_true_2, 'black', ...
    T, y_b_1, 'blue', T, y_b_2, 'blue')