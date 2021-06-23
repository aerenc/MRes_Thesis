% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

clc;
clear;

global x_jmt beta alpha xi_jmt Total J M T index index2 index3 inits 

rand('seed',2344)    % reset uniform random number generator
randn('seed',2344)   % reset normal random number generator

%% SIMULATION SETTING - Generating Data

% Generating the data for logit case first:

alpha   = 2;                    % true alpha
beta    = [-6,3,0.5,0.5]';      % true betas
gamma   = [5,0.5,0.5,0.5]';     % true gamma
theta   = [alpha; beta; gamma]; % true theta

I       = 200;                  % Total number of consumers in each market
J       = 20;                   % Initial total good + outside share in each market
M       = 20;                   % Total markets
T       = 10;                   % Total time period
Total   = J*M*T;                % Total goods present a priori

index   = [J:J:Total]';             % indexing each good in every market M*T combination
indexx  = [J-1:J-1:Total-M*T]';     % indexing that will be used in calculating outside good share in this setting

index2             = zeros(M*T,2);  % indexing the inside goods in each market in matrix form
index2(1,1)        = 1;
index2(1,2)        = J-1;
for i = 2 : M*T
    index2(i,1)    = index(i-1,1)+1;       
    index2(i,2)    = index(i,1)-1; 
end
index3             = zeros(M*T,2);  % indexing optimizer prices
index3(1,1)        = 1   ;
index3(1,2)        = 19;
for i = 2 : M*T
   index3(i,1)     = index3(i-1,2)+1;
   index3(i,2)     = index3(i,1) + 18;
end


x_jmt_1 = ones(Total,1);
x_jmt_1(index,:) = 0;
x_jmt_2 = unifrnd(1,5,Total,1);
x_jmt_2(index,:) = 0;
x_jmt_3 = unifrnd(0,1,Total,1);
x_jmt_3(index,:) = 0;
x_jmt_4 = unifrnd(-2,3,Total,1);
x_jmt_4(index,:) = 0;
x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3 x_jmt_4]; % Generated observed product characteristics

w_jmt_1 = ones(Total,1);
w_jmt_1(index,:) = 0;
w_jmt_2 = unifrnd(0,1,Total,1);
w_jmt_2(index,:) = 0;
w_jmt_3 = unifrnd(-2,3,Total,1);
w_jmt_3(index,:) = 0;
w_jmt_4 = unifrnd(0.1,0.9,Total,1);
w_jmt_4(index,:) = 0;
w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3 w_jmt_4]; % Generated observed cost shifters


% Now, it's time to construct our unobserved product characteristics and
% unobserved costs. We will introduce correlation across markets, across
% time periods, and across the same good:

x       = normrnd(0,1,J-1,1);
m       = normrnd(0,0.2,M,1);
t       = normrnd(0,0.4,T,1);


xi_jmt  = zeros(J,M,T);

for   j = 1:J-1
    for market = 1:M
       for time = 1:T
            xi_jmt(j,market,time) = x(j,1) + m(market,1) + t(time,1);
       end
    end
end

xi_jmt             = xi_jmt(:);                 % unobserved product characteristics vector
xi_jmt(index,:)    = 0;

omega_jmt          = normrnd(0,1,Total,1);  % unobserved costs
omega_jmt(index,:) =0;

mc_jmt             = w_jmt * gamma + omega_jmt; % getting marginal costs
mc__jmt            = zeros(Total-M*T,1);

for i = 1: M*T
    
mc__jmt(index3(i,1):index3(i,2))  = mc_jmt(index2(i,1):index2(i,2),1);
end



inits        = ones(Total-M*T,1);

solveforprices     = @(ppp) (ppp -(1./(alpha.*(1 - s_jmt(ppp)))) - mc__jmt);
 % ppp -> 3800*1; s_jmt(ppp) ->4000*1: boyutlar bu handle'da tutmuyor. 
 
 
 

%opts        = optimset('Algorithm ','interior-point ','Display ','iter ');
%opts        = optimset('Display','off');
opts         = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000); %


tic
p_jmt        = fsolve(solveforprices,inits,opts); % this gives us prices
toc



s_jmt          = s_jmt(p_jmt);          % gives us market shares

% array_s_jmt  = reshape(s_jmt,J,M,T);  % market share in array form
% array_p_jmt  = reshape(p_jmt,J,M,T);  % prices in array form
% array_xi_jmt = reshape(xi_jmt,J,M,T); % unobs. prod. chars. in array form


cumsum_s_jmt         = cumsum(s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 

s0_jmt               = 1 - marketwisesum_s_jmt;











