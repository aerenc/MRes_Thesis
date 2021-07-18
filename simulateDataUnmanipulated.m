function [x_jmt, w_jmt, xi_jmt, omega_jmt, mc_jmt] = simulateDataUnmanipulated(J,M,T,Total,gamma)

global x m t extra_noise
% excluding the outside good <=> Total - M*T

x_jmt_1 = unifrnd(0,1,Total-M*T,1);
x_jmt_2 = unifrnd(0,1,Total-M*T,1);
x_jmt_3 = unifrnd(0,1,Total-M*T,1);
x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3];                                       % Generated observed product characteristics

w_jmt_1 = unifrnd(0,1,Total-M*T,1);
w_jmt_2 = unifrnd(0,1,Total-M*T,1);
w_jmt_3 = unifrnd(0,1,Total-M*T,1);
w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3];                                       % Generated cost shifters

% Now, it's time to construct our unobserved product characteristics and
% unobserved costs. We will introduce correlation across markets, across
% time periods, and across the same good:

x       = normrnd(1,3,J-1,1);                                            % Brand-specific components
m       = normrnd(0.5,2.4,M,1);                                            % Market-fixed effects
t       = normrnd(0.5,1.8,T,1);                                            % Year-fixed effects

xi_jmt  = zeros(J-1,M,T);

extra_noise = normrnd(0,1,J-1,M,T);                                        % Extra random shock

for   j = 1:J-1
  for market = 1:M
    for time = 1:T
        xi_jmt(j,market,time) = x(j,1) + m(market,1) + t(time,1) + extra_noise(j,market,time);
    end
  end
end

xi_jmt               = xi_jmt(:);                                          % Unobserved product characteristics vector
omega_jmt            = normrnd(0,1,Total-M*T,1);                           % Unobserved costs (drawn from normal distribution with no particular correlation structure)
mc_jmt               = [ones(Total-M*T,1) w_jmt] * gamma + omega_jmt;      % Marginal cost specification

end