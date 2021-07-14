%% RM THESIS IN EOR - 2021

%% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

%%

clc;
clear;

global  beta alpha Total I J M T   

randn('seed',10);               % Reset normal random number generator
rand('seed',10);                % Reset uniform random number generator
%% SIMULATION SETTING - LOGIT

% Generating the data for logit case first, we will extend the required data for the RC logit:

alpha   = 2;                    % True alpha
beta    = [4,3,0.5,0.5]';       % True betas
gamma   = [5,1,0.5,0.5]';       % True gamma
sigma   = [0.5,0.8]';           % Heterogeneity term for RC logit
theta   = [alpha; beta; gamma]; % True theta

I       = 200;                  % Total number of consumers in each market (for estimation part)
J       = 6;                    % Initial total good + outside share in each market
M       = 30;                   % Total markets
T       = 10;                   % Total time period
Total   = J*M*T;                % Total goods present a priori, INCLUDING outside good

%% CONSTRUCTING INDEXES:

% We construct some indexes for the simulated data: straight out indexes
% correspond to the unmanipulated data and m_ indexes correspond to the
% manipulated data, also constructing the dropping index which is defined
% as in indexes function:

[index ,indexx ,index2 ,index3] = indexes(J,M,T,Total); 

%% SIMULATING THE BASE DATA:

% Herein, I simulate the baseline data:

[base_x_jmt, base_w_jmt, base_xi_jmt, base_omega_jmt, base_mc_jmt] = simulateDataUnmanipulated(J,M,T,Total,gamma);


%% SOLVING FOR SHARES & PRICES FOR UNMANIPULATED DATA:

inits                = ones(Total-M*T,1);                                  % Getting initial values for getting eq. prices


solveforprices       = @(ppp) (ppp -(1./(alpha.*(1 - shares(ppp,base_x_jmt,base_xi_jmt,index,index2,index3)))) - base_mc_jmt);
opts                 = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000);

tic
base_p_jmt           = fsolve(solveforprices,inits,opts);                  % This gives us equilibrium prices for logit-based simulation
toc

base_s_jmt           = shares(base_p_jmt,base_x_jmt,base_xi_jmt,index,index2,index3);% this gives us inner market shares

cumsum_s_jmt         = cumsum(base_s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
base_s0_mt           = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo

clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt marketwisesum_s_jmt

%% UNMANIPULATED DATA IN 3D-ARRAY ((J-1)*M*T) FORM AND MATRIX FORM (reshaped: dimensions:(J-1)*(M*T)):

[base_array_s_jmt,base_array_p_jmt,base_array_xi_jmt,base_array_omega_jmt,base_array_mc_jmt,base_array_x_jmt,base_array_w_jmt,   ...
base_reshaped_x1_jmt,base_reshaped_x2_jmt,base_reshaped_x3_jmt,base_reshaped_w1_jmt,base_reshaped_w2_jmt,base_reshaped_w3_jmt] = ...
reshapeDataUnmanipulated(J,M,T,base_s_jmt,base_p_jmt,base_xi_jmt,base_omega_jmt,base_mc_jmt,base_x_jmt,base_w_jmt);

% %% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA & OLS ESTIMATION OF FIRST STAGE IV FOR MANIPULATED DATA & 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA
% 
% [m_OLSestimates,m_xi_jmt_error,m_Twoslsestimatessecondstage] = OLS_IV_estimation_Manipulated(Total,M,T,adj,m_index3,m_s_jmt,m_s0_mt,m_x_jmt,m_p_jmt,...
%  m_w_jmt);
% 
%% SIMULATION SETTING - RC LOGIT

% We will introduce heterogeneity on constant term and price:

v        = randn(2,I);                                                     % Draws for share integrals during estimation (for constant + price)
sigma    = [0.5,0.8]';                                                     % Heterogeneity term for RC logit (for constant + price)
theta_RC = [alpha;beta;sigma;gamma];                                       % True theta


inits_RC             = ones(Total-M*T,1);                                  % Initial values for price-share optimization

solveforprices_RC    = @(pp) (pp -(1./(alpha.*(1 - shares_RC(pp,base_x_jmt,base_xi_jmt,sigma,index,index2,index3)))) - base_mc_jmt);
opts                 = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000);

tic
base_p_jmt_RC        = fsolve(solveforprices_RC,inits,opts);               % This gives us equilibrium prices for RC-based simulation
toc

base_s_jmt_RC        = shares_RC(base_p_jmt_RC,base_x_jmt,base_xi_jmt,sigma,index,index2,index3);     % This gives us inner market shares

cumsum_s_jmt         = cumsum(base_s_jmt_RC);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
base_s0_mt_RC             = 1 - marketwisesum_s_jmt;                       % Getting outside shares for each M*T combo

clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt marketwisesum_s_jmt


