% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER
clc;
clear;
global x_jmt beta alpha xi_jmt Total J M T index index2 index3 inits adj m_x_jmt m_xi_jmt m_index m_index2 m_index3 adj
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

% Constructing the BLP instruments: Note that we will have two types of
% instruments: own product characteristics and characteristics sum of rivals' characteristics 
% across markets and time periods

cumsum_z1_jmt     = cumsum(x_jmt);
cumsum_z1         = cumsum_z1_jmt(J-1,:);
totalsum_z1       = diff(cumsum_z1_jmt(indexx,:)); 
totalsum_z1       = [cumsum_z1;totalsum_z1]; %(M*T)*4 -> contains the sum of product characteristics for every M*T combination

for  i = 1: M*T
     rival_jmt(index3(i,1):index3(i,2),:) = totalsum_z1(i,:) - x_jmt(index3(i,1):index3(i,2),:);   % sum of rivals' observed product characteristics 
end

Z_blp = [ones(Total-M*T,1) x_jmt rival_jmt x_jmt.^2 rival_jmt.^2]; %BLP instruments -> contains own & rivals prod. characteristics and their squares.


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


%% DATA MANIPULATION

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



indexxx = [1:1:J-1]';                        % indexing from 1 to 19 (i.e. max. number of brands for all M*T)

adj = 5;
willbedroppedelementindex = zeros(adj,M*T);  % these elements will be dropped

manipulatingindex1 = randi([1 4],1,M*T);
manipulatingindex2 = randi([5 8],1,M*T);
manipulatingindex3 = randi([9 12],1,M*T);
manipulatingindex4 = randi([13 16],1,M*T);
manipulatingindex5 = randi([17 19],1,M*T);

willbedroppedelementindex(1,:) =  manipulatingindex1;
willbedroppedelementindex(2,:) =  manipulatingindex2;
willbedroppedelementindex(3,:) =  manipulatingindex3;
willbedroppedelementindex(4,:) =  manipulatingindex4;
willbedroppedelementindex(5,:) =  manipulatingindex5;  

reshaped_x1_jmt = reshape(x_jmt(:,1),J-1,M*T); 
reshaped_x2_jmt = reshape(x_jmt(:,2),J-1,M*T); 
reshaped_x3_jmt = reshape(x_jmt(:,3),J-1,M*T); 

reshaped_w1_jmt = reshape(w_jmt(:,1),J-1,M*T); 
reshaped_w2_jmt = reshape(w_jmt(:,2),J-1,M*T); 
reshaped_w3_jmt = reshape(w_jmt(:,3),J-1,M*T); 

reshaped_xi_jmt    = reshape(xi_jmt,J-1,M*T); 
reshaped_omega_jmt = reshape(omega_jmt,J-1,M*T);

for  i = 1 : M*T
    reshaped_x1_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_x2_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_x3_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_w1_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_w2_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_w3_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_xi_jmt(willbedroppedelementindex(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_omega_jmt(willbedroppedelementindex(:,i),i) = 0;
end

%%
% m_ stands for manipulated:

m_x1_jmt = reshaped_x1_jmt(:);
m_x1_jmt = m_x1_jmt(m_x1_jmt~=0);
m_x2_jmt = reshaped_x2_jmt(:);
m_x2_jmt = reshaped_x1_jmt(m_x2_jmt~=0);
m_x3_jmt = reshaped_x3_jmt(:);
m_x3_jmt = m_x3_jmt(m_x3_jmt~=0);

m_x_jmt = [m_x1_jmt m_x2_jmt m_x3_jmt] ;

m_w1_jmt = reshaped_w1_jmt(:);
m_w1_jmt = m_w1_jmt(m_w1_jmt~=0);
m_w2_jmt = reshaped_w2_jmt(:);
m_w2_jmt = m_w2_jmt(m_w2_jmt~=0);
m_w3_jmt = reshaped_w3_jmt(:);
m_w3_jmt = m_w3_jmt(m_w3_jmt~=0);

m_w_jmt = [m_w1_jmt m_w2_jmt m_w3_jmt];

m_xi_jmt = reshaped_xi_jmt(:);
m_xi_jmt = m_xi_jmt(m_xi_jmt~=0);

m_omega_jmt = reshaped_omega_jmt(:);
m_omega_jmt = m_omega_jmt(m_omega_jmt~=0);

m_mc_jmt = [ones(Total-M*T-adj*M*T,1) m_w_jmt] * gamma + m_omega_jmt; % getting marginal costs for manipulated data



m_index   = [J-adj:J-adj:Total-adj*M*T]';             % indexing each good in every market M*T combination for manipulated case
m_indexx  = [J-1-adj:J-1-adj:Total-M*T-adj*M*T]';     

m_index2             = zeros(M*T,2);             % indexing the inside goods in each market in matrix form for manipulated case
m_index2(1,1)        = 1;
m_index2(1,2)        = J-1-adj;
for i = 2 : M*T
    m_index2(i,1)    = m_index(i-1,1)+1;       
    m_index2(i,2)    = m_index(i,1)-1; 
end
m_index3             = zeros(M*T,2);  % indexing optimizer prices
m_index3(1,1)        = 1;
m_index3(1,2)        = 14;
for i = 2 : M*T
   m_index3(i,1)     = m_index3(i-1,2)+1;
   m_index3(i,2)     = m_index3(i,1) + 13;
end

m_inits              = ones(Total-M*T-adj*M*T,1);
m_solveforprices     = @(pp) (pp -(1./(alpha.*(1 - m_s_jmt(pp)))) - m_mc_jmt);
%opts        = optimset('Algorithm ','interior-point ','Display ','iter ');
%opts        = optimset('Display','off');
opts                 = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000); %
tic
m_p_jmt              = fsolve(m_solveforprices,m_inits,opts); % this gives us prices
toc

m_s_jmt                = m_s_jmt(m_p_jmt);          % gives us inner market shares
cumsum_m_s_jmt         = cumsum(m_s_jmt);
cumsum_m_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_m_s_jmt  = diff(cumsum_m_s_jmt(m_indexx,:)); 
marketwisesum_m_s_jmt  = [cumsum_m_s_jmt_1;marketwisesum_m_s_jmt]; 
m_s0_jmt               = 1 - marketwisesum_m_s_jmt;





