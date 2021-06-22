% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

clc;
clear;

global x_jmt beta alpha xi_jmt Total J

rand('seed',2344)    % reset uniform random number generator
randn('seed',2344)   % reset normal random number generator

%% SIMULATION SETTING - Generating Data

% Generating the data for logit case first:

% alpha = 2;                 % true alpha
% beta  = [3,3,0.5,0.5]';         % true betas
% gamma = [5,0.5,0.5,0.5]';         % true gamma
% theta = [alpha; beta; gamma]; % true theta
% 
% I     = 200;                  % Total number of consumers in each market
% J     = 15;                   % Initial total good in each market; some observations will be removed
% M     = 20;                   % Total markets
% T     = 10;                   % Total time period
% 
% Total = J*M*T;                % Total goods present a priori
% 
% x_jmt_1 = ones(Total,1);
% x_jmt_2 = unifrnd(15,200,Total,1);
% x_jmt_3 = unifrnd(0,1,Total,1);
% x_jmt_4 = unifrnd(1,50,Total,1);
% x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3 x_jmt_4]; % Generated observed product characteristics
% 
% w_jmt_1 = ones(Total,1);
% w_jmt_2 = unifrnd(10,100,Total,1);
% w_jmt_3 = unifrnd(0,1,Total,1);
% w_jmt_4 = unifrnd(1,10,Total,1);
% w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3 w_jmt_4]; % Generated observed cost shifters

alpha = -2.5;                 % true alpha
beta  = [3,2,5,1.5]';         % true betas
gamma = [5,1,3,0.5]';         % true gamma
theta = [alpha; beta; gamma]; % true theta

I     = 200;                  % Total number of consumers in each market
J     = 15;                   % Initial total good in each market; some observations will be removed
M     = 20;                   % Total markets
T     = 10;                   % Total time period

Total = J*M*T;                % Total goods present a priori

x_jmt_1 = ones(Total,1);
x_jmt_2 = unifrnd(3,7,Total,1);
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
omega_jmt = normrnd(0,1,Total,1);      % unobserved costs

mc_jmt    = w_jmt * gamma + omega_jmt; % getting marginal costs


solveforprices = @(ppp) (ppp - (1./(alpha.*(1 - s_jmt(ppp)))) - mc_jmt);


%opts   = optimset('Display','off');
 opts    = optimset('Display','iter','TolCon',1E-6,'TolFun',1E-6,'TolX',1E-10,'MaxFunEvals',1000000000); %

tic
p_jmt  = fsolve(solveforprices,ones(Total,1),opts); % this gives us prices
toc
% p_jmt = fminunc(solveforprices,ones(Total,1),opts);

s_jmt = s_jmt(p_jmt);    % gives us market shares



