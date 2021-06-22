% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

clc;
clear;

rand('seed',2344)    % reset uniform random number generator
randn('seed',2344)   % reset normal random number generator

%% SIMULATION SETTING

% Generating the data for logit case first:

alpha = -2.5;                 % true alpha
beta  = [4,5,2,1.5]';       % true betas
gamma = [3,1,1,0.5]';         % true gamma
theta = [alpha; beta; gamma]; % true theta

I     = 200;                  % Total number of consumers in each market
J     = 15;                   % Initial total good in each market; some observations will be removed
M     = 20;                   % Total markets
T     = 10;                   % Total time period

Total = J*M*T;                % Total goods present a priori

x_jmt_1 = ones(Total,1);
x_jmt_2 = unifrnd(4,6,Total,1);
x_jmt_3 = unifrnd(-1,5,Total,1);
x_jmt_4 = unifrnd(0.5,2.5,Total,1);
x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3 x_jmt_4]; % Generated observed product characteristics

w_jmt_1 = ones(Total,1);
w_jmt_2 = unifrnd(0,2,Total,1);
w_jmt_3 = unifrnd(0.5,1.5,Total,1);
w_jmt_4 = unifrnd(0.1,0.9,Total,1);
w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3 w_jmt_4]; % Generated observed cost shifters


% Now, it's time to construct our unobserved product characteristics and
% unobserved costs. We will introduce correlation across markets, across
% time periods, and across the same good:

x = normrnd(0,1,J,1);
m = normrnd(0,0.2,M,1);
t = normrnd(0,0.4,T,1);


xi_jmt = zeros(J,M,T);

for j = 1:J
    for market = 1:M
        for time = 1:T
            xi_jmt(j,market,time) = x(j,1) + m(market,1) + t(time,1);
        end
    end
end

xi_jmt    = xi_jmt(:);                 % unobserved product characteristics vector
omega_jmt = normrnd(0,0.5,Total,1);    % unobserved costs

mc_jmt    = w_jmt * gamma + omega_jmt; % getting marginal costs

p_jmt = rand(Total,1);                % prices (will be calculated)
% 
% tic
% p_jmt = SolveforPrice_logit(x_jmt,beta,alpha,xi_jmt,w_jmt,gamma,omega_jmt,Total);
% toc

delta_mt          = x_jmt * beta + alpha .* p_jmt + xi_jmt; 
nom               = exp(delta_mt);
cumsum_delta_jmt  = cumsum(nom);
cumsum_1          = cumsum_delta_jmt(15,:);
index             = [15:15:3000]';
totalsum          = diff(cumsum_delta_jmt(index,:)); 
totalsum          = [cumsum_1;totalsum];  % mt total share of each j \in J_{mt}
Totalsum          = repelem(totalsum,15); % stretching it out to 3000*1 (since we have 15 products in each market a priori)
denom             = ones(Total,1) + Totalsum;
s_jmt             = nom./denom;


solveforprices = @(p_jmt) (log(p_jmt - (ones(Total,1)./(alpha.*(ones(Total,1) - s_jmt)))) - w_jmt * gamma - omega_jmt);

%solveforprices = @(p) (log(p - (1./(alpha.*(1 - s_mt./(1+ repelem([cumsum(exp(x_jmt * beta + p .* alpha + xi_jmt))(15,:); diff(cumsum(exp(x_jmt*beta+p.*alpha+xi_jmt))([15:15:3000]',:))],15))))))) - w_jmt * gamma - omega_jmt);



%opts      = optimset('Display','off');
opts    = optimset('Display','iter','TolCon',1E-16,'TolFun',1E-16,'TolX',1E-16);

tic
prices  = fsolve(solveforprices,rand(Total,1),opts);
toc








