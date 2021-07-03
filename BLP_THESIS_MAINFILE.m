% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

clc;
clear;

global  beta alpha Total J M T   

rng('default');

%% SIMULATION SETTING - LOGIT

% Generating the data for logit case first:

alpha   = 2;                    % true alpha
beta    = [5,3,0.5,0.5]';       % true betas
gamma   = [5,1,0.5,0.5]';       % true gamma
theta   = [alpha; beta; gamma]; % true theta

I       = 200;                  % Total number of consumers in each market
J       = 20;                   % Initial total good + outside share in each market
M       = 20;                   % Total markets
T       = 10;                   % Total time period
Total   = J*M*T;                % Total goods present a priori, INCLUDING outside good
adj     = 5;                    % 5 goods will be removed from each market & time period randomly as in manipulated case

%% INDEXES:

[index ,indexx ,index2 ,index3 ,m_index ,m_indexx ,m_index2 ,m_index3, m_drop_index] = indexes(J,M,T,Total,adj);

%% SIMULATING THE DATA:

[x_jmt ,w_jmt ,xi_jmt ,omega_jmt ,mc_jmt] = simulateDataUnmanipulated(J,M,T,Total,gamma);


%% SOLVING FOR SHARES & PRICES:

inits                = ones(Total-M*T,1);
solveforprices       = @(ppp) (ppp -(1./(alpha.*(1 - shares(ppp,x_jmt,xi_jmt,index,index2,index3,0)))) - mc_jmt);
opts                 = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000);

tic
p_jmt                = fsolve(solveforprices,inits,opts);                    % this gives us prices
toc

s_jmt                = shares(p_jmt,x_jmt,xi_jmt,index,index2,index3,0);     % this gives us inner market shares
cumsum_s_jmt         = cumsum(s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt                = 1 - marketwisesum_s_jmt;                              % getting outside shares for each M*T combo

%% UNMANIPULATED DATA IN 3D-ARRAY ((J-1)*M*T) FORM AND MATRIX FORM (dimensions:(J-1)*(M*T)):

[array_s_jmt,array_p_jmt,array_xi_jmt,array_omega_jmt,array_mc_jmt,array_x_jmt,array_w_jmt,reshaped_x1_jmt,reshaped_x2_jmt,reshaped_x3_jmt,...
 reshaped_w1_jmt,reshaped_w2_jmt,reshaped_w3_jmt] = reshapeDataUnmanipulated(J,M,T,s_jmt,p_jmt,xi_jmt,omega_jmt,mc_jmt,x_jmt,w_jmt);


%% MANIPULATION BEGINS -> GETTING MANIPULATED SIMULATED DATA

[m_x_jmt,m_w_jmt,m_xi_jmt,m_omega_jmt,m_mc_jmt] = simulateDataManipulated(Total,J,M,T,adj,x_jmt,w_jmt,xi_jmt,omega_jmt,m_drop_index,gamma);


%% SOLVING FOR SHARES & PRICES FOR MANIPULATED CASE:

m_inits               = ones(Total-M*T-adj*M*T,1);
m_solveforprices      = @(pp) (pp -(1./(alpha.*(1 - shares(pp,m_x_jmt,m_xi_jmt,m_index,m_index2,m_index3,adj)))) - m_mc_jmt);

tic
m_p_jmt               = fsolve(m_solveforprices,m_inits,opts);                            % this gives us prices
toc

m_s_jmt               = shares(m_p_jmt,m_x_jmt,m_xi_jmt,m_index,m_index2,m_index3,adj);   % this gives us inner market shares
cumsum_m_s_jmt        = cumsum(m_s_jmt);
cumsum_m_s_jmt_1      = cumsum_m_s_jmt(J-adj-1,:);
marketwisesum_m_s_jmt = diff(cumsum_m_s_jmt(m_indexx,:)); 
marketwisesum_m_s_jmt = [cumsum_m_s_jmt_1;marketwisesum_m_s_jmt]; 
m_s0_mt               = 1 - marketwisesum_m_s_jmt;              % this gives us market shares of outside good in each market and time period

%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA & OLS ESTIMATION OF FIRST STAGE IV FOR MANIPULATED DATA & 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA

[m_OLSestimates,m_xi_jmt_error,m_Twoslsestimatessecondstage] = OLS_IV_estimation_Manipulated(Total,M,T,adj,m_index3,m_s_jmt,m_s0_mt,m_x_jmt,m_p_jmt,...
 m_w_jmt);

%% INTRODUCTION OF EXISTING PRODUCT IN LOGIT CASE











