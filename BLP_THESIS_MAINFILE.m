%% RM THESIS IN EOR - 2021

%% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC-LOGIT SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

%%

clc;
clear;

global  beta alpha Total I J M T x m t extra_noise

rng('default');

%% 1 - SIMULATION:

% In the first part of this code, I will generate the required data:

%% SIMULATION SETTING - LOGIT

% Generating the data for logit case first, we will extend the required data for the RC logit:

alpha         = 2;                    % True alpha
beta          = [5,3,0.5,0.5]';       % True betas (B0 ; B1 ; B2 ; B3)
gamma         = [5,1,0.5,0.5]';       % True gamma
theta_logit   = [alpha; beta; gamma]; % True theta (but w/o heterogeneity (i.e. sigma) parameters)

I             = 200;                  % Total number of consumers in each market (will be used in estimation part)
J             = 11;                   % Initial total good + outside share in each market
M             = 30;                   % Total markets
T             = 10;                   % Total time period
Total         = J*M*T;                % Total goods present "a priori" i.e. "real" data + first good in first market at first time period,
                                      % INCLUDING outside goods for each M*T market-time combo

%% CONSTRUCTING INDEXES:

% We construct some indexes for the simulated data: straight out indexes
% correspond to the unmanipulated data and m_ indexes correspond to the
% manipulated data, also constructing the dropping index which is defined
% as in indexes function:

[base_index, base_indexx, base_index2 ,base_index3, index, indexx, index2, index3] = indexes(J,M,T,Total); 

%% SIMULATING THE BASE DATA:

% Herein, I simulate the baseline data:

[base_x_jmt, base_w_jmt, base_xi_jmt, base_omega_jmt, base_mc_jmt] = simulateDataUnmanipulated(J,M,T,Total,gamma);
x;
m;
t;
extra_noise;

%% SOLVING FOR SHARES & PRICES FOR UNMANIPULATED DATA BY LOGIT SPECIFICATION:

base_inits                 = ones(Total-M*T,1);                            % Getting initial values for getting eq. prices

base_solveforprices_logit  = @(ppp) (ppp -(1./(alpha.*(1 - base_shares(ppp,base_x_jmt,base_xi_jmt,base_index,base_index2,base_index3)))) - base_mc_jmt);
opts                 = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000);

tic
base_p_jmt_logit     = fsolve(base_solveforprices_logit,base_inits,opts);  % This gives us equilibrium prices for logit-based simulation
toc

base_s_jmt_logit     = base_shares(base_p_jmt_logit,base_x_jmt,base_xi_jmt,base_index,base_index2,base_index3);
                                                                           % This gives us inner market shares                                                                         
cumsum_s_jmt         = cumsum(base_s_jmt_logit);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(base_indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
base_s0_mt_logit     = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt marketwisesum_s_jmt

%% UNMANIPULATED DATA IN 3D-ARRAY ((J-1)*M*T) FORM AND MATRIX FORM (reshaped: dimensions:(J-1)*(M*T)):

[base_array_s_jmt,base_array_p_jmt,base_array_xi_jmt,base_array_omega_jmt,base_array_mc_jmt,base_array_x_jmt,base_array_w_jmt,   ...
base_reshaped_x1_jmt,base_reshaped_x2_jmt,base_reshaped_x3_jmt,base_reshaped_w1_jmt,base_reshaped_w2_jmt,base_reshaped_w3_jmt] = ...
reshapeDataUnmanipulated(J,M,T,base_s_jmt_logit,base_p_jmt_logit,base_xi_jmt,base_omega_jmt,base_mc_jmt,base_x_jmt,base_w_jmt);
 
%% SIMULATION SETTING - RC-LOGIT:

% We will introduce heterogeneity on constant term and price:

v        = randn(2,I);                                                     % Draws for share integrals during estimation (for constant + price)
sigma    = [0.5,0.8]';                                                     % Heterogeneity term for RC logit (for constant + price)
theta    = [alpha;beta;sigma;gamma];                                       % True theta

%% SOLVING FOR SHARES & PRICES FOR UNMANIPULATED DATA BY RC-LOGIT SPECIFICATION:

base_solveforprices = @(pp) (pp -(1./(alpha.*(1 - base_shares_RC(pp,base_x_jmt,base_xi_jmt,sigma,base_index,base_index2,base_index3)))) - base_mc_jmt);

tic
base_p_jmt             = fsolve(base_solveforprices,base_inits,opts);      % This gives us equilibrium prices for RC-based simulation
toc

base_s_jmt             = base_shares_RC(base_p_jmt,base_x_jmt,base_xi_jmt,sigma,base_index,base_index2,base_index3);     
                                                                           % This gives us inner market shares
cumsum_s_jmt           = cumsum(base_s_jmt);
cumsum_s_jmt_1         = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt    = diff(cumsum_s_jmt(base_indexx,:)); 
marketwisesum_s_jmt    = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
base_s0_mt             = 1 - marketwisesum_s_jmt;                          % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt 

%% MANIPULATION: GETTING THE "REAL" DATA:

% Now, we will drop the first good from the first market in the first time
% period and resimulate the data for logit and RC-logit cases:

[x_jmt, w_jmt, xi_jmt, omega_jmt, mc_jmt]  = simulateDataManipulated(base_x_jmt,base_w_jmt,base_xi_jmt,base_omega_jmt,base_mc_jmt);

%% SOLVING FOR SHARES & PRICES FOR "REAL" DATA USING LOGIT SPECIFICATION:

inits               = ones(Total-M*T-1,1);                                 % Getting initial values for getting eq. prices

solveforprices_logit= @(pr) (pr -(1./(alpha.*(1 - shares(pr,x_jmt,xi_jmt,index,index2,index3)))) - mc_jmt);

tic
p_jmt_logit         = fsolve(solveforprices_logit,inits,opts);             % This gives us equilibrium prices for logit-based simulation
toc

s_jmt_logit         = shares(p_jmt_logit,x_jmt,xi_jmt,index,index2,index3);% This gives us inner market shares

cumsum_s_jmt         = cumsum(s_jmt_logit);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt_logit          = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt

%% SOLVING FOR SHARES & PRICES FOR "REAL" DATA USING RC-LOGIT SPECIFICATION:

solveforprices       = @(pr) (pr -(1./(alpha.*(1 - shares_RC(pr,x_jmt,xi_jmt,sigma,index,index2,index3)))) - mc_jmt);

tic
p_jmt                = fsolve(solveforprices,inits,opts);                  % This gives us equilibrium prices for logit-based simulation
toc

s_jmt                = shares_RC(p_jmt,x_jmt,xi_jmt,sigma,index,index2,index3);   
                                                                           % This gives us inner market shares
cumsum_s_jmt         = cumsum(s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt                = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt 

%% 2 - ESTIMATION:

% Until this point, I generated the whole required baseline data, did some
% simulation to arrive at equilibrium "data" prices and market shares. 

%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR "REAL" DATA & OLS ESTIMATION OF FIRST STAGE IV FOR "REAL" DATA & 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR "REAL" DATA:

one = 1;

[OLSestimates_logit,xi_jmt_error_logit,Twoslsestimatessecondstage_logit]   = Estimation_OLS_IV(Total,M,T,one,index3,s_jmt_logit,s0_mt_logit,x_jmt,p_jmt_logit,w_jmt); % GIVING COMPLEX ESTIMS
[OLSestimates,xi_jmt_error,Twoslsestimatessecondstage]                     = Estimation_OLS_IV(Total,M,T,one,index3,s_jmt,s0_mt,x_jmt,p_jmt,w_jmt);


%% "BLP" (RC-LOGIT) ESTIMATION OF THE MODEL:

REALDATA  = [x_jmt w_jmt p_jmt s_jmt];                                     % THIS GIVES ALL REQUIRED DATA FOR BLP ESTIMATION








