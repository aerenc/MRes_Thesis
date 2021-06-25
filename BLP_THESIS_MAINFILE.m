% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER
clc;
clear;
global x_jmt beta alpha xi_jmt Total J M T index index2 index3 inits 
rng('default');
%% SIMULATION SETTING - Generating Data
% Generating the data for logit case first:

alpha   = 2;                    % true alpha
beta    = [5,3,0.5,0.5]';       % true betas
gamma   = [5,1,0.5,0.5]';       % true gamma
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

x_jmt_1 = unifrnd(0,1,Total-M*T,1);
x_jmt_2 = unifrnd(0,1,Total-M*T,1);
x_jmt_3 = unifrnd(0,1,Total-M*T,1);
x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3]; % Generated observed product characteristics

w_jmt_1 = unifrnd(0,1,Total-M*T,1);
w_jmt_2 = unifrnd(0,1,Total-M*T,1);
w_jmt_3 = unifrnd(0,1,Total-M*T,1);
w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3]; % Generated observed cost shifters

% Now, it's time to construct our unobserved product characteristics and
% unobserved costs. We will introduce correlation across markets, across
% time periods, and across the same good:

x       = normrnd(0,1,J-1,1);
m       = normrnd(0,0.2,M,1);
t       = normrnd(0,0.4,T,1);

xi_jmt  = zeros(J-1,M,T);

for   j = 1:J-1
    for market = 1:M
       for time = 1:T
            xi_jmt(j,market,time) = x(j,1) + m(market,1) + t(time,1);
       end
    end
end

xi_jmt             = xi_jmt(:);                 % unobserved product characteristics vector
omega_jmt          = normrnd(0,1,Total-M*T,1);      % unobserved costs
mc_jmt             = [ones(Total-M*T,1) w_jmt] * gamma + omega_jmt; % getting marginal costs
inits              = ones(Total-M*T,1);
solveforprices     = @(ppp) (ppp -(1./(alpha.*(1 - s_jmt(ppp)))) - mc_jmt);
%opts        = optimset('Algorithm ','interior-point ','Display ','iter ');
%opts        = optimset('Display','off');
opts         = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000); %
tic
p_jmt        = fsolve(solveforprices,inits,opts); % this gives us prices
toc

s_jmt                = s_jmt(p_jmt);          % gives us inner market shares
cumsum_s_jmt         = cumsum(s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_jmt               = 1 - marketwisesum_s_jmt;

% Constructing the instruments: Note that we will have two types of
% instruments: own product characteristics and characteristics sum of rivals' characteristics 
% across markets and time periods



cumsum_z1_jmt     = cumsum(x_jmt);
cumsum_z1         = cumsum_z1_jmt(J-1,:);
totalsum_z1       = diff(cumsum_z1_jmt(indexx,:)); 
totalsum_z1       = [cumsum_z1;totalsum_z1]; %(M*T)*4 -> contains the sum of product characteristics for every M*T combination

for  i = 1: M*T
     rival_jmt(index3(i,1):index3(i,2),:) = totalsum_z1(i,:) - x_jmt(index3(i,1):index3(i,2),:);   % sum of rivals' observed product characteristics 
end

Z = [ones(Total-M*T,1) x_jmt rival_jmt x_jmt.^2 rival_jmt.^2]; %instruments -> contains own & rivals prod. characteristics and their squares.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now we need to construct log(S_{jm}/S_{0m}) dependent variable:

dependentlogratio = zeros(Total-M*T,1);

for i = 1 : M*T
    dependentlogratio(index3(i,1):index3(i,2),:) = log(s_jmt(index3(i,1):index3(i,2),:))-log(s0_jmt(i,1));
end

X            = [ones(Total-M*T,1) x_jmt p_jmt];              % stacking regressors together
OLSestimates = (X'*X)\X'*dependentlogratio; % classical OLS formula giving us OLS estimates on aggregate Logit model


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF FIRST STAGE IV %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, to be able to remove the endogeneity bias in our model, we regress
% price additionally on cost shifters and collect the corresponding estimates for prices

Z = [ones(Total-M*T,1) w_jmt x_jmt];                                 
price_hat = Z*inv(Z'*Z)*Z'*p_jmt;               % part of price that has no correlation with error term
Twoslsestimatesfirststage = inv(Z'*Z)*Z'*p_jmt; % First stage estimates of coefficients



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since we have the estimates for prices, we can now estimate our model by
% also estimating the second stage of IV regression:

X_hat = [ones(Total-M*T,1) x_jmt price_hat];
Twoslsestimatessecondstage = inv(X_hat' * X_hat)*X_hat'*dependentlogratio;


%% 24 june 2021 - gereken her sey generate edildi; simdi data elimination yapmak gerekiyor

% putting all required data in array form:

array_s_jmt     = reshape(s_jmt,J-1,M,T);         % market share in array form
array_p_jmt     = reshape(p_jmt,J-1,M,T);         % prices in array form
array_xi_jmt    = reshape(xi_jmt,J-1,M,T);       % unobs. prod. chars. in array form
array_omega_jmt = reshape(omega_jmt,J-1,M,T);
array_mc_jmt    = reshape(mc_jmt,J-1,M,T);

array_x_jmt     = zeros(J-1,M,T,size(x_jmt,2));
x11_jmt = x_jmt(:,1);
x12_jmt = x_jmt(:,2);
x13_jmt = x_jmt(:,3);
array_x_jmt(:,:,:,1) = reshape(x11_jmt,J-1,M,T); % unobs. prod. chars. in array form
array_x_jmt(:,:,:,2) = reshape(x12_jmt,J-1,M,T); % unobs. prod. chars. in array form
array_x_jmt(:,:,:,3) = reshape(x13_jmt,J-1,M,T); % unobs. prod. chars. in array form

array_w_jmt = zeros(J-1,M,T,size(w_jmt,2));
w11_jmt = w_jmt(:,1);
w12_jmt = w_jmt(:,2);
w13_jmt = w_jmt(:,3);
array_w_jmt(:,:,:,1) = reshape(w11_jmt,J-1,M,T); % unobs. prod. chars. in array form
array_w_jmt(:,:,:,2) = reshape(w12_jmt,J-1,M,T); % unobs. prod. chars. in array form
array_w_jmt(:,:,:,3) = reshape(w13_jmt,J-1,M,T); % unobs. prod. chars. in array form

% Now comes the data manipulation. We will drop some observations across
% markets and time periods.

% Intuition: Delete the good in each market (excluding market M=20)
% corresponding the the INDEX of that market.
% E.g. : delete good 1 from m=1 for every t
       % delete good 2 from m=2 for every t, iterate this process


indexxx = [1:1:J-1]';        % indexing from 1 to 19 (i.e. max. number of brands for all M*T)

willbedroppedelementindex = zeros(5*M,T);  % these elements will be dropped

manipulatingindex1 = randi([1 4],M,T);
manipulatingindex2 = randi([5 8],M,T);
manipulatingindex3 = randi([9 12],M,T);
manipulatingindex4 = randi([13 16],M,T);
manipulatingindex5 = randi([17 19],M,T);


willbedroppedelementindex(1:5:5*M-4,:) =  manipulatingindex1;
willbedroppedelementindex(2:5:5*M-3,:) =  manipulatingindex2;
willbedroppedelementindex(3:5:5*M-2,:) =  manipulatingindex3;
willbedroppedelementindex(4:5:5*M-1,:) =  manipulatingindex4;
willbedroppedelementindex(5:5:5*M,:)   =  manipulatingindex5;  












