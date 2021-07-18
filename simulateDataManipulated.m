function [x_jmt, w_jmt, xi_jmt, omega_jmt, mc_jmt]  = simulateDataManipulated(base_x_jmt,base_w_jmt,base_xi_jmt, base_omega_jmt, base_mc_jmt)

% Excluding the outside good                  <=> Total - M*T
% Excluding the outside good + j,m,t = 1,1,1  <=> Total - 1 - M*T

% => Since we already simulated these before, we will just drop the first
% data point from each:

x_jmt_1 = base_x_jmt(2:end,1);
x_jmt_2 = base_x_jmt(2:end,2);
x_jmt_3 = base_x_jmt(2:end,3);
x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3];                                       % Generated observed product characteristics

w_jmt_1 = base_w_jmt(2:end,1);
w_jmt_2 = base_w_jmt(2:end,2);
w_jmt_3 = base_w_jmt(2:end,3);
w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3];                                       % Generated cost shifters

xi_jmt  = base_xi_jmt(2:end,:);

omega_jmt            = base_omega_jmt(2:end,:);                            % Unobserved costs (drawn from normal distribution with no particular correlation structure)
mc_jmt               = base_mc_jmt(2:end,:);                               % Marginal cost specification

end




