%% RM THESIS IN EOR - 2021

%% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC-LOGIT SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

%%

clc;
clear;

global  beta alpha Total I J M T x m t extra_noise constant share IV price indexxx W Ktheta IDmkt ns A z tol_inner Kbeta v one Indicator error error_J iterativedelta Kmc zero marginalcost g gg zeta heterogeneity base_heterogeneity unobs_cost alpha_i simulatedshare RCderivative price_updated share_updated error_J_updated gg_updated marginalcost_updated zeta_updated unobs_cost_updated alpha_i_updated simulatedshare_updated RCderivative_updated beta_est alpha_est sigma_est

rng('default');                                                            % Setting seed for the whole thesis

%% 1 - SIMULATION:

% In the first part of this code, I will generate the required data:

%% SIMULATION SETTING - LOGIT

% Generating the data for logit case first, we will extend the required data for the RC logit:

alpha         = 2;                    % True alpha
beta          = [5,3,0.5,0.5]';       % True betas (B0 ; B1 ; B2 ; B3)
gamma         = [5,0.5,0.5,0.5]';     % True gamma
theta_logit   = [beta;-alpha;gamma];  % True theta (but w/o heterogeneity (i.e. sigma) parameters)

I             = 30;                   % Total number of consumers in each market (will be used in estimation part)
J             = 6;                    % Initial total good + outside share in each market
M             = 20;                   % Total markets
T             = 10;                   % Total time period
Total         = J*M*T;                % Total goods present "a priori" i.e. "real" data + first good in first market at first time period,
                                      % INCLUDING outside goods for each M*T market-time combo
one           = 1;
%% CONSTRUCTING INDEXES:

% We construct some indexes for the simulated data: straight out indexes
% correspond to the unmanipulated data and m_ indexes correspond to the
% manipulated data, also constructing the dropping index which is defined
% as in indexes function:

[base_index, base_indexx, base_index2 ,base_index3, index, indexx, index2, index3, indexxx, IDmkt] = indexes(J,M,T,Total); 

%% SIMULATING THE BASE DATA:

% Herein, I simulate the baseline data:

[base_x_jmt, base_w_jmt, base_xi_jmt, base_omega_jmt, base_mc_jmt] = simulateDataUnmanipulated(J,M,T,Total,gamma);
x;
m;
t;
extra_noise;

%% SOLVING FOR SHARES & PRICES FOR UNMANIPULATED DATA BY LOGIT SPECIFICATION:

base_inits                 = zeros(Total-M*T,1);                           % Getting initial values for getting eq. prices

base_solveforprices_logit  = @(ppp) (ppp -(1./(alpha.*(1 - base_shares(ppp,base_x_jmt,base_xi_jmt,base_index,base_index2,base_index3)))) - base_mc_jmt);
opt       = optimset('Display','iter','TolCon',1E-6,'TolFun',1E-10,'TolX',1E-10);

tic
base_p_jmt_logit     = fsolve(base_solveforprices_logit,base_inits,opt);  % This gives us equilibrium prices for logit-based simulation
toc

base_s_jmt_logit     = base_shares(base_p_jmt_logit,base_x_jmt,base_xi_jmt,base_index,base_index2,base_index3);
                                                                           % This gives us inner market shares                                                                         
cumsum_s_jmt         = cumsum(base_s_jmt_logit);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(base_indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
base_s0_mt_logit     = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt marketwisesum_s_jmt

%% SIMULATION SETTING - RC-LOGIT:

% We will introduce heterogeneity on constant term and price:

sigmaa   = [0.5]';                                                     % Heterogeneity term for RC logit (for price)
theta    = [alpha;beta;sigmaa;gamma];                                      % True theta
base_heterogeneity = randn(1,I);

%% SOLVING FOR SHARES & PRICES FOR UNMANIPULATED DATA BY RC-LOGIT SPECIFICATION:

base_solveforprices = @(pp) (pp -(1./(alpha.*(1 - base_shares_RC(pp,base_x_jmt,base_xi_jmt,sigmaa,base_index,base_index2,base_index3)))) - base_mc_jmt);

tic
base_p_jmt             = fsolve(base_solveforprices,base_inits,opt);      % This gives us equilibrium prices for RC-based simulation
toc

base_s_jmt             = base_shares_RC(base_p_jmt,base_x_jmt,base_xi_jmt,sigmaa,base_index,base_index2,base_index3);     
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

inits                = ones(Total-M*T-1,1);                                 % Getting initial values for getting eq. prices

solveforprices_logit = @(pr) (pr -(1./(alpha.*(1 - shares(pr,x_jmt,xi_jmt,index,index2,index3)))) - mc_jmt);

tic
p_jmt_logit         = fsolve(solveforprices_logit,inits,opt);             % This gives us equilibrium prices for logit-based simulation
toc

s_jmt_logit         = shares(p_jmt_logit,x_jmt,xi_jmt,index,index2,index3);% This gives us inner market shares

cumsum_s_jmt         = cumsum(s_jmt_logit);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt_logit          = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt

%% SOLVING FOR SHARES & PRICES FOR "REAL" DATA USING RC-LOGIT SPECIFICATION:

%heterogeneity        = randn(1,I);
heterogeneity = [0.472787993096300,-1.16798166954839,2.31864830172844,-0.997857136400435,-0.704987439045016,0.907164087520108,0.175537518078233,0.0682443228063713,0.00591625073786253,-1.33957733079912,-1.28545090125272,-1.67652301869667,0.755868016253087,-0.895768868332429,-1.73464079444528,0.476550199886526,-0.253969587070142,-1.30315683387076,0.290260392687052,0.439663451333156,0.503393235574352,-0.610059855304868,-0.353936356966121,-0.432099688896280,0.164478755253172,0.510496995046600,1.64276118009046,-0.615126067924010,-0.628671339470476,-0.176557908312171];

solveforprices       = @(pr) (pr -(1./(alpha.*(1 - shares_RC(pr,x_jmt,xi_jmt,sigmaa,index,index2,index3)))) - mc_jmt);

tic
p_jmt                = fsolve(solveforprices,inits,opt);                   % This gives us equilibrium prices for RC logit-based simulation
toc

s_jmt                = shares_RC(p_jmt,x_jmt,xi_jmt,sigmaa,index,index2,index3);   
                                                                           % This gives us inner market shares
cumsum_s_jmt         = cumsum(s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt                = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt 

%% 2 - ESTIMATION:

% Until this point, I generated the whole required baseline data, did some
% simulation to arrive at equilibrium "data" prices and market shares. Now it's
% time to estimate the constructed model:

%% A) OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR "REAL" DATA & OLS ESTIMATION OF FIRST STAGE IV FOR "REAL" DATA & 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR "REAL" DATA:

%[OLSestimates,xi_jmt_error,Twoslsestimatessecondstage] = Estimation_OLS_IV(Total,M,T,one,index3,s_jmt,s0_mt,x_jmt,p_jmt,w_jmt);
                                                                           % This gives both OLS + IV estimates for "real" data
                                                                                                                                                      
[OLSestimates,xi_jmt_error_logit,Twoslsestimatessecondstage,xi_jmt_error_twosls,se_ols,se_twosls] =...
    Estimation_OLS_IV(Total,M,T,one,index3,s_jmt,s0_mt,x_jmt,p_jmt,w_jmt);                                                                           
                                                                           
%% B) "BLP" (RC-LOGIT) ESTIMATION OF THE MODEL:

%% B1) DEMAND ESTIMATION:

REALDATA  = [x_jmt w_jmt p_jmt s_jmt];                                     % THIS GIVES ALL REQUIRED DATA FOR BLP ESTIMATION

A     = REALDATA(:,1:3);
z     = REALDATA(:,4:6);
price = REALDATA(:,7);
share = REALDATA(:,8);

constant  = ones(Total-M*T-1,1);                                           % Creating constant for regression   
tol_inner = 1.e-14;                                                        % Tolerance for inner loop (NFXP)
ns        = 50;
Kbeta     = 2+size(A,2);                                                   % =5(constant&price&prod. characteristics) - # of parameters in mean utility
Ktheta    = 1;                                                             % =2(price&constant)          - # of parameters with random coefficient
%v         = randn(Ktheta,ns);                                              % Draws for share integrals during estimation
v         = [1.16125413303262,-0.398468092854212,0.928355195202805,-0.735626270177556,-0.855994624611678,1.81311089390539,0.860349737621205,0.567824302200727,0.00664181484955421,-0.599786854214078,1.17523115009483,1.13125585851945,-1.76191503658499,-1.15183604026249,-0.368293013835837,1.48413803681476,0.391096577027212,1.30765291645093,-0.577178565951493,0.477100296010847,-0.621523889494184,-1.87514401626840,0.307786196845486,0.0390454996297821,-0.761705086537191,-0.549914774549615,-0.779440665402206,0.520463114879734,-1.04080086146455,0.591864830585222,-0.291608736237927,-0.236180217045201,-0.876779018252948,-0.836816341912153,-0.672945800076526,-0.919695072404640,-0.706421954056373,-0.642052446655979,-1.51973772137986,-1.74539386354885,-3.12890661533416,0.545687916708685,-0.974727752231523,0.221010185690586,-0.394221229858108,0.837841052967890,0.379404226606850,1.47959824333850,0.238174938734941,0.538982839646241];                                               % Draws for share integrals during estimation

IV        = [constant A z A.^2 z.^2];                                      % Instruments
%IV       = [constant A z];                                 

nIV       = size(IV,2);                                                    % # of instrumental variables
W         = (IV'*IV)\eye(nIV);                                             % Starting GMM weighting matrix
init      = rand(Kbeta+Ktheta,1);                                          % Starting values for all parameters; serves to simulate them uniformly

x_L       = [-Inf.*ones(Kbeta,1);zeros(Ktheta,1)];                         % Lower bounds is zero for standard deviations of random coefficients--
                                                                           % (for theta_2 parameters)                                                                         
x_U       = [Inf.*ones(Kbeta,1);Inf.*ones(Ktheta,1)];                      % Upper bounds for standard deviations of random coefficients--
                                                                           % (for theta_2 parameters)

iterativedelta = zeros(Total-M*T-1,1);                                     % Creating a global variable for delta which will be updated each time 
                                                                           % in contraction mapping until a fixed point is attained

Indicator      = zeros(M*T,Total-M*T-1);                                   % Matrix whose rows consist of markets and columns consist of all products: 
                                                                           % It basically assigns value "1" for products present in a given market, 
                                                                           % and "0" otherwise. 
for k = 1:M*T
    Indicator(k,index3(k,1):index3(k,2)) = 1;
end
%% % FOR J=6;M=20;T=10;I=50:

tic
[X,fval_rep,exitflag,output,lambda,grad,hessian] = fmincon(@BLP_GMM_demand,init,[],[],[],[],x_L,x_U,[],opt);
toc

demandmoment = g;

xi_D         = error;

theta1D = X(1:Kbeta,1);                                                    % Estimated Mean Tastes

theta2D = X(Kbeta+1:Kbeta+Ktheta,1);                                       % Estimated Deviations from Mean Taste

Grad1 = zeros(size(demandmoment,1),Kbeta+Ktheta);
for i = 1:Kbeta+Ktheta
    Grad1(:,i) = gradient(demandmoment, X(i,1));                           % Getting gradients of estimates to be used in variance matrix
end
varianceest_d = pinv(Grad1'*W*Grad1)*Grad1'*W*(demandmoment*demandmoment')*W*Grad1*pinv(Grad1'*W*Grad1); % using GMM variance matrix formulae
sterrors_d    = sqrt(diag(varianceest_d));                                 % Getting standard errors of estimated parameters

%% B2) JOINT ESTIMATION:

Kmc     = 1+size(z,2);                                                     % =4 (constant & cost shifters) # of parameters in MC function
zero    = zeros(size(W,1),size(W,1));                                      % Getting matrix filled with 0s in dimension of W to be used in joint moment minimization

%init_J  = rand(Kbeta+Ktheta+Kmc,1);                                        % Initial values for this part
init_J  = [0.427772303103763;0.567792073333150;0.613038459346904;0.806993254654286;0.880783416271628;0.159978307447974;0.0891769254299106;0.313134600651602;0.607352744347650;0.378738421986569];                                        % Initial values for this part

x_L_J     = [-Inf.*ones(Kbeta,1);zeros(Ktheta,1);-Inf.*ones(Kmc,1)];       % Lower bounds  
x_U_J     = Inf.*ones(Kbeta+Ktheta+Kmc,1);                                 % Upper bounds 

tic
[X_J,fval_rep2,exitflag2,output2,lambda2] = fmincon(@BLP_GMM_joint,init_J,[],[],[],[],x_L_J,x_U_J,[],opt);
toc

xi_J           = error_J;

alpha_RC = alpha_i;

s_ijmt  = simulatedshare;

ds_dp = RCderivative;

demandmoment_J = gg;

supplymoment_J = zeta;

omega_J = unobs_cost;

mc_J           = marginalcost;

theta1J = X_J(1:Kbeta,1);                        % Estimated Mean Tastes

theta2J = X_J(Kbeta+1:Kbeta+Ktheta,1);           % Estimated Deviations from Mean Taste

gammaJ  = X_J(Kbeta+Ktheta+1:Kbeta+Ktheta+Kmc);  % MC estimates

theta_J = [theta1J; theta2J];

Grad2 = zeros(size(demandmoment_J,1),Kbeta+Ktheta);
for i = 1:Kbeta+Ktheta
    Grad2(:,i) = gradient(demandmoment_J, X_J(i,1));                       % Getting gradients of demand-side estimates to be used in variance matrix
end
varianceestQ2D  = pinv(Grad2'*W*Grad2)*Grad2'*W*(demandmoment_J*demandmoment_J')*W*Grad2*pinv(Grad2'*W*Grad2); 
                                                                           % Using GMM variance matrix formulae for demand-side
sterrorsQ2D     = sqrt(diag(varianceestQ2D));                              % Getting standard errors of demand-side

Grad3 = zeros(size(supplymoment_J,1),Kmc);
for i = 1:Kmc
    Grad3(:,i) = gradient(supplymoment_J, gammaJ(i,1));                    % Getting gradients of supply-side estimates to be used in variance matrix
end
varianceestQ2S = pinv(Grad3'*W*Grad3)*Grad3'*W*(supplymoment_J*supplymoment_J')*W*Grad3*pinv(Grad3'*W*Grad3);  
                                                                           % Using GMM variance matrix formulae for supply-side
sterrorsQ2S     = sqrt(diag(varianceestQ2S));                              % Getting standard errors of supply-side
sterrors_j      = [sterrorsQ2D; sterrorsQ2S];                              % Stacking demand-side & supply-side errors together

%% SOLVING FOR SHARES & PRICES FOR "REAL" DATA USING RC-LOGIT SPECIFICATION: (UPDATED)
% HERE, I USE DS/DP I GET FROM FIRST BLP JOINT ESTIMATION TO REGENERATE
% PRICES:

solveforprices_updated       = @(pr) (pr - (shares_RC(pr,x_jmt,xi_jmt,sigmaa,index,index2,index3)./(-ds_dp))  - mc_jmt);

tic
p_jmt_updated                = fsolve(solveforprices_updated,inits,opt);                   % This gives us equilibrium prices for RC logit-based simulation
toc

s_jmt_updated                = shares_RC(p_jmt_updated,x_jmt,xi_jmt,sigmaa,index,index2,index3);   
                                                                           % This gives us inner market shares
cumsum_s_jmt         = cumsum(s_jmt_updated);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt_updated        = 1 - marketwisesum_s_jmt;                            % Getting outside shares for each M*T combo
clear cumsum_s_jmt cumsum_s_jmt_1 marketwisesum_s_jmt 


%% B2) JOINT ESTIMATION "UPDATED":

price_updated = p_jmt_updated;
share_updated = s_jmt_updated;

tic
[X_J_updated,fval_rep2_updated,exitflag2_updated,output2_updated,lambda2_updated] = fmincon(@BLP_GMM_joint_updated,init_J,[],[],[],[],x_L_J,x_U_J,[],opt);
toc

xi_J_updated           = error_J_updated;

alpha_RC_updated = alpha_i_updated;

s_ijmt_updated  = simulatedshare_updated;

ds_dp_updated = RCderivative_updated;

demandmoment_J_updated = gg_updated;

supplymoment_J_updated = zeta_updated;

omega_J_updated = unobs_cost_updated;

mc_J_updated           = marginalcost_updated;

theta1J_updated = X_J_updated(1:Kbeta,1);                        % Estimated Mean Tastes

theta2J_updated = X_J_updated(Kbeta+1:Kbeta+Ktheta,1);           % Estimated Deviations from Mean Taste

gammaJ_updated  = X_J_updated(Kbeta+Ktheta+1:Kbeta+Ktheta+Kmc);  % MC estimates

theta_J_updated = [theta1J_updated; theta2J_updated];

Grad2_updated = zeros(size(demandmoment_J_updated,1),Kbeta+Ktheta);
for i = 1:Kbeta+Ktheta
    Grad2_updated(:,i) = gradient(demandmoment_J_updated, X_J_updated(i,1));                       % Getting gradients of demand-side estimates to be used in variance matrix
end
varianceestQ2D_updated  = pinv(Grad2_updated'*W*Grad2_updated)*Grad2_updated'*W*(demandmoment_J_updated*demandmoment_J_updated')*W*Grad2_updated*pinv(Grad2_updated'*W*Grad2_updated); 
                                                                           % Using GMM variance matrix formulae for demand-side
sterrorsQ2D_updated     = sqrt(diag(varianceestQ2D_updated));                              % Getting standard errors of demand-side

Grad3_updated = zeros(size(supplymoment_J_updated,1),Kmc);
for i = 1:Kmc
    Grad3_updated(:,i) = gradient(supplymoment_J_updated, gammaJ_updated(i,1));                    % Getting gradients of supply-side estimates to be used in variance matrix
end
varianceestQ2S_updated = pinv(Grad3_updated'*W*Grad3_updated)*Grad3_updated'*W*(supplymoment_J_updated*supplymoment_J_updated')*W*Grad3_updated*pinv(Grad3_updated'*W*Grad3_updated);  
                                                                           % Using GMM variance matrix formulae for supply-side
sterrorsQ2S_updated     = sqrt(diag(varianceestQ2S_updated));                              % Getting standard errors of supply-side
sterrors_j_updated      = [sterrorsQ2D_updated; sterrorsQ2S_updated];                              % Stacking demand-side & supply-side errors together


%%

x_111 = base_x_jmt(1,:);
w_111 = base_w_jmt(1,:);

xxx            = zeros(Total-M*T-one,1);         % Brand dummies
xxx(1:J-2,:)   = [2 3 4 5]';
xxx(J-1:end,:) = repmat([1 2 3 4 5]',M*T-1,1);

mmm                  = zeros(Total-M*T-one,1);   % Market dummies
mmm(1:J-2,:)         = [1 1 1 1]';
mmm(J-1:(J-1)*M-1,:) = repelem([2:1:20]',(J-1),1);
mmm((J-1)*M:end,:)   = repmat(repelem([1:1:M]',(J-1),1),(T-1),1);

ttt                  = zeros(Total-M*T-one,1);   % Time dummies
ttt(1:(J-1)*M-1,:)   = 1 ;
ttt((J-1)*M:end,:)   = repelem ([2:1:T]',(J-1)*M,1);

X_xi      = [xxx mmm ttt];                                                 % Stacking regressors together
estims_xi = (X_xi'*X_xi)\X_xi'*xi_J;                                       % Classical OLS formula giving us regression estimates for xi regression that researcher does
wololoo   = xi_J - X_xi * estims_xi ;                                      % This gives estimated unobserved product characteristics
nnn=size(X_xi,1);
kkk=size(X_xi,2);
cov_xi=(wololoo'*wololoo/(nnn-kkk))*inv(X_xi'*X_xi);  
se_xi=sqrt(diag(cov_xi));                                               % This identifies \sigma_x, \sigma_m and \sigma_t

first_index = [J-1:J-1:(J-1)*(M*T-1)]';
second_index = [J-1:J-1:(J-1)*(M-1)]';


first_xi = xi_J(first_index,:);
second_xi = xi_J(second_index,:);
third_xi = normrnd(estims_xi(1,:),se_xi(1,:),100,1) + normrnd(estims_xi(2,:),se_xi(2,:),100,1) + normrnd(estims_xi(3,:),se_xi(3,:),100,1);

first_xi_mean = mean(first_xi);
second_xi_mean = mean(second_xi);
third_xi_mean = mean(third_xi);

beta_est  = X_J(1:4,:);
alpha_est = X_J(5,:);
sigma_est = X_J(6,:);
gamma_est = X_J(7:end,:);
mc_est    = [[1 w_111]*gamma_est+base_omega_jmt(1,:); mc_J(1:J-2)];  

inits_xi = ones(J-1,1);

for k = 1: size(first_xi,1)
    
solveforprices_first_xi    = @(pppp) (pppp -(1./((-alpha_est).*(1 - shares_RC_counterfactual(pppp,base_x_jmt,xi_J,first_xi(k,:))))) - mc_est);
    
tic
p_jmt_first_xi(k,:)                = fsolve(solveforprices_first_xi,inits_xi,opt);             
toc

s_jmt_first_xi(k,:) = shares_RC_counterfactual(p_jmt_first_xi(k,:)',base_x_jmt,xi_J_updated,first_xi(k,:));

end


for k = 1: size(second_xi,1)
    
solveforprices_second_xi    = @(pppp) (pppp -(1./((-alpha_est).*(1 - shares_RC_counterfactual(pppp,base_x_jmt,xi_J,second_xi(k,:))))) - mc_est);
    
tic
p_jmt_second_xi(k,:)                = fsolve(solveforprices_second_xi,inits_xi,opt);             
toc

s_jmt_second_xi(k,:) = shares_RC_counterfactual(p_jmt_second_xi(k,:)',base_x_jmt,xi_J,second_xi(k,:));

end


for k = 1: size(third_xi,1)
    
solveforprices_third_xi    = @(ppppp) (ppppp -(1./((-alpha_est).*(1 - shares_RC_counterfactual(ppppp,base_x_jmt,xi_J,third_xi(k,:))))) - mc_est);
    
tic
p_jmt_third_xi(k,:)                = fsolve(solveforprices_third_xi,inits_xi,opt);             
toc

s_jmt_third_xi(k,:) = shares_RC_counterfactual(p_jmt_third_xi(k,:)',base_x_jmt,xi_J,third_xi(k,:));

end














































































































% x_111 = base_x_jmt(1,:);
% w_111 = base_w_jmt(1,:);
% 
% xxx            = zeros(Total-M*T-one,1);         % Brand dummies
% xxx(1:J-2,:)   = [2 3 4 5]';
% xxx(J-1:end,:) = repmat([1 2 3 4 5]',M*T-1,1);
% 
% mmm                  = zeros(Total-M*T-one,1);   % Market dummies
% mmm(1:J-2,:)         = [1 1 1 1]';
% mmm(J-1:(J-1)*M-1,:) = repelem([2:1:20]',(J-1),1);
% mmm((J-1)*M:end,:)   = repmat(repelem([1:1:M]',(J-1),1),(T-1),1);
% 
% ttt                  = zeros(Total-M*T-one,1);   % Time dummies
% ttt(1:(J-1)*M-1,:)   = 1 ;
% ttt((J-1)*M:end,:)   = repelem ([2:1:T]',(J-1)*M,1);
% 
% X_xi      = [xxx mmm ttt];                                                 % Stacking regressors together
% estims_xi = (X_xi'*X_xi)\X_xi'*xi_J;                                       % Classical OLS formula giving us regression estimates for xi regression that researcher does
% wololoo   = xi_J - X_xi * estims_xi ;                                      % This gives estimated unobserved product characteristics
% nnn=size(X_xi,1);
% kkk=size(X_xi,2);
% cov_xi=(wololoo'*wololoo/(nnn-kkk))*inv(X_xi'*X_xi);  
% se_xi=sqrt(diag(cov_xi));                                               % This identifies \sigma_x, \sigma_m and \sigma_t
% 
% first_index = [J-1:J-1:(J-1)*(M*T-1)]';
% second_index = [J-1:J-1:(J-1)*(M-1)]';
% 
% 
% first_xi = xi_J(first_index,:);
% second_xi = xi_J(second_index,:);
% third_xi = normrnd(estims_xi(1,:),se_xi(1,:),100,1) + normrnd(estims_xi(2,:),se_xi(2,:),100,1) + normrnd(estims_xi(3,:),se_xi(3,:),100,1);
% 
% first_xi_mean = mean(first_xi);
% second_xi_mean = mean(second_xi);
% third_xi_mean = mean(third_xi);
% 
% beta_est  = X_J(1:4,:);
% alpha_est = X_J(5,:);
% sigma_est = X_J(6,:);
% gamma_est = X_J(7:end,:);
% mc_est    = [[1 w_111]*gamma_est+base_omega_jmt(1,:); mc_J(1:J-2)];  
% 
% inits_xi = 20.*ones(J-1,1);
% 
% for k = 1: size(first_xi,1)
%     
% solveforprices_first_xi    = @(pppp) (pppp -(1./((-alpha_est).*(1 - shares_RC_counterfactual(pppp,base_x_jmt,xi_J,first_xi(k,:))))) - mc_est);
%     
% tic
% p_jmt_first_xi(k,:)                = fsolve(solveforprices_first_xi,inits_xi,opt);             
% toc
% 
% s_jmt_first_xi(k,:) = shares_RC_counterfactual(p_jmt_first_xi(k,:)',base_x_jmt,xi_J,first_xi(k,:));
% 
% end
% 
% for k = 1: size(second_xi,1)
%     
% solveforprices_second_xi    = @(pppp) (pppp -(1./((-alpha_est).*(1 - shares_RC_counterfactual(pppp,base_x_jmt,xi_J,second_xi(k,:))))) - mc_est);
%     
% tic
% p_jmt_second_xi(k,:)                = fsolve(solveforprices_second_xi,inits_xi,opt);             
% toc
% 
% s_jmt_second_xi(k,:) = shares_RC_counterfactual(p_jmt_second_xi(k,:)',base_x_jmt,xi_J,second_xi(k,:));
% 
% end
% 
% for k = 1: size(third_xi,1)
%     
% solveforprices_third_xi    = @(ppppp) (ppppp -(1./((-alpha_est).*(1 - shares_RC_counterfactual(ppppp,base_x_jmt,xi_J,third_xi(k,:))))) - mc_est);
%     
% tic
% p_jmt_third_xi(k,:)                = fsolve(solveforprices_third_xi,inits_xi,opt);             
% toc
% 
% s_jmt_third_xi(k,:) = shares_RC_counterfactual(p_jmt_third_xi(k,:)',base_x_jmt,xi_J,third_xi(k,:));
% 
% end

























