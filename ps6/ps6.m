%PS6

t_0 = 0;
y_0 = 2;
t_y = t_0:250;
t_z = t_0:5:250;

Prob2_single_state_truth_gen;
hold on;
plot(ytrue);
scatter(tmeas,z);

figure
plot(utrue)
hold on

Prob2_dual_state_truth_gen
[y5,t5] = ...
    Prob2_dual_state_lin_model( ...
    y0true,alpha,utrue,Nsteps,dt,t_beg,w_variance);
plot(y5(2,:))
alpda(4) = 3;
[y3,t3] = ...
    Prob2_dual_state_lin_model( ...
    y0true,alpha,utrue,Nsteps,dt,t_beg,w_variance);
plot(y3(2,:))
alpda(4) = 7;
[y7,t7] = ...
    Prob2_dual_state_lin_model( ...
    y0true,alpha,utrue,Nsteps,dt,t_beg,w_variance);
plot(y7(2,:))
alpda(4) = 1;
[y1,t1] = ...
    Prob2_dual_state_lin_model( ...
    y0true,alpha,utrue,Nsteps,dt,t_beg,w_variance);
plot(y1(2,:))
alpda(4) = 9;
[y9,t9] = ...
    Prob2_dual_state_lin_model( ...
    y0true,alpha,utrue,Nsteps,dt,t_beg,w_variance);
plot(y9(2,:))