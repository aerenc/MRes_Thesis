%% RM THESIS IN EOR - 2021

%% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC-LOGIT SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

%%

clc;
clear;

global  beta alpha Total I J M T x m t extra_noise constant share IV price indexxx W Ktheta IDmkt ns A z tol_inner Kbeta v one Indicator error error_J iterativedelta Kmc zero marginalcost g gg zeta heterogeneity unobs_cost alpha_i beta_est alpha_est sigma_est base_x_jmt xi_J

%rng('default');                                                           % Setting seed for the whole thesis
rand('seed',88);
randn('seed',44);
%% 1 - SIMULATION:

% In the first part of this code, I will generate the required data:

%% SIMULATION SETTING - LOGIT

% Generating the data for logit case first, we will extend the required data for the RC logit:

alpha         = 2;                                                         % True alpha
beta          = [5,3,0.5,0.5]';                                            % True betas (B0 ; B1 ; B2 ; B3)
gamma         = [5,0.5,0.5,0.5]';                                          % True gamma

I             = 50;                                                        % Total number of consumers in each market (will be used in estimation part)
J             = 6;                                                         % Initial total good + outside share in each market
M             = 20;                                                        % Total markets
T             = 10;                                                        % Total time period
Total         = J*M*T;                                                     % Total goods present "a priori" i.e. "real" data + first good in first market at first time period,
                                                                           % INCLUDING outside goods for each M*T market-time combo
one           = 1;

opt      = optimset('Display','iter','TolCon',1E-6,'TolFun',1E-10,'TolX',1E-10);
sigmaa   = 0.5;                                                            % Heterogeneity term for RC logit (for price)
theta    = [beta;-alpha;sigmaa;gamma];                                     % True theta
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

%% MANIPULATION: GETTING THE "REAL" DATA:

% Now, we will drop the first good from the first market in the first time
% period and resimulate the data for logit and RC-logit cases:

[x_jmt, w_jmt, xi_jmt, omega_jmt, mc_jmt]  = simulateDataManipulated(base_x_jmt,base_w_jmt,base_xi_jmt,base_omega_jmt,base_mc_jmt);

%% SOLVING FOR SHARES & PRICES FOR "REAL" DATA USING RC-LOGIT SPECIFICATION:

inits                = ones(Total-M*T-1,1);                                % Getting initial values for getting eq. prices
heterogeneity        = randn(1,I);

solveforprices       = @(pr) (pr - shares_RC(pr,x_jmt,xi_jmt,sigmaa,index,index2,index3) - mc_jmt);

tic
p_jmt                = fsolve(solveforprices,inits,opt);                   % This gives us equilibrium prices for RC logit-based simulation
toc

s_jmt                = sharescalculator_RC(p_jmt,x_jmt,xi_jmt,sigmaa,index,index2,index3);   
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

[OLSestimates,xi_jmt_error_logit,Twoslsestimatessecondstage,xi_jmt_error_twosls,se_ols,se_twosls] =...
    Estimation_OLS_IV(Total,M,T,one,index3,s_jmt,s0_mt,x_jmt,p_jmt,w_jmt);    
                                                                           % This gives both OLS + IV estimates for "real" data
                                                                                                                                                                                    
%% B) "BLP" (RC-LOGIT) ESTIMATION OF THE MODEL:

%% B1) DEMAND ESTIMATION:

REALDATA  = [x_jmt w_jmt p_jmt s_jmt];                                     % THIS GIVES ALL REQUIRED DATA FOR BLP ESTIMATION

A     = REALDATA(:,1:3);
z     = REALDATA(:,4:6);
price = REALDATA(:,7);
share = REALDATA(:,8);

constant  = ones(Total-M*T-1,1);                                           % Creating constant for regression   
tol_inner = 1.e-14;                                                        % Tolerance for inner loop (NFXP)
ns        = 100;
Kbeta     = 2+size(A,2);                                                   % =5(constant&price&prod. characteristics) - # of parameters in mean utility
Ktheta    = 1;                                                             % =2(price&constant)          - # of parameters with random coefficient
v         = randn(Ktheta,ns);                                              % Draws for share integrals during estimation
IV        = [constant A z A.^2 z.^2];                                      % Instruments

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

clear Grad1 varianceest_d demandmoment xi_D
%% B2) JOINT ESTIMATION:

Kmc     = 1+size(z,2);                                                     % =4 (constant & cost shifters) # of parameters in MC function
zero    = zeros(size(W,1),size(W,1));                                      % Getting matrix filled with 0s in dimension of W to be used in joint moment minimization

init_J  = rand(Kbeta+Ktheta+Kmc,1);                                        % Initial values for this part

x_L_J     = [-Inf.*ones(Kbeta,1);zeros(Ktheta,1);-Inf.*ones(Kmc,1)];       % Lower bounds  
x_U_J     = Inf.*ones(Kbeta+Ktheta+Kmc,1);                                 % Upper bounds 

tic
[X_J,fval_rep2,exitflag2,output2,lambda2] = fmincon(@BLP_GMM_joint,init_J,[],[],[],[],x_L_J,x_U_J,[],opt);
toc

xi_J     = error_J;

alpha_RC = alpha_i;

demandmoment_J = gg;

supplymoment_J = zeta;

omega_J        = unobs_cost;

mc_J           = marginalcost;

theta1J        = X_J(1:Kbeta,1);                                           % Estimated Mean Tastes

theta2J        = X_J(Kbeta+1:Kbeta+Ktheta,1);                              % Estimated Deviations from Mean Taste

gammaJ         = X_J(Kbeta+Ktheta+1:Kbeta+Ktheta+Kmc);                     % MC estimates

theta_J        = [theta1J; theta2J];

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

clear Grad2 Grad3 varianceestQ2D sterrorsQ2D varianceestQ2S sterrorsQ2S s_ijmt ds_dp demandmoment_J supplymoment_J
%% 3 - COUNTERFACTUAL ANALYSIS:

% Here, I insert the "missing" first good in M=1&T=1 to that market itself
% to check new prices, shares, profits, and CS:

x_111 = base_x_jmt(1,:);
w_111 = base_w_jmt(1,:);

xxx            = zeros(Total-M*T-one,1);                                   % Brand dummies
xxx(1:J-2,:)   = [2 3 4 5]';
xxx(J-1:end,:) = repmat([1 2 3 4 5]',M*T-1,1);

mmm                  = zeros(Total-M*T-one,1);                             % Market dummies
mmm(1:J-2,:)         = [1 1 1 1]';
mmm(J-1:(J-1)*M-1,:) = repelem([2:1:20]',(J-1),1);
mmm((J-1)*M:end,:)   = repmat(repelem([1:1:M]',(J-1),1),(T-1),1);

ttt                  = zeros(Total-M*T-one,1);                             % Time dummies
ttt(1:(J-1)*M-1,:)   = 1 ;
ttt((J-1)*M:end,:)   = repelem ([2:1:T]',(J-1)*M,1);

X_xi      = [xxx mmm ttt];                                                 % Stacking regressors together
estims_xi = (X_xi'*X_xi)\X_xi'*xi_J;                                       % Classical OLS formula giving us regression estimates for xi regression that researcher does
wololoo   = xi_J - X_xi * estims_xi ;                                      
nnn=size(X_xi,1);
kkk=size(X_xi,2);
cov_xi=(wololoo'*wololoo/(nnn-kkk))*inv(X_xi'*X_xi);  
se_xi=sqrt(diag(cov_xi));                                                  % This identifies \sigma_x, \sigma_m and \sigma_t

first_index  = [J-1:J-1:(J-1)*(M*T-1)]';
second_index = [J-1:J-1:(J-1)*(M-1)]';


first_xi    = xi_J(first_index,:);
second_xi   = xi_J(second_index,:);


% Now I use the regression to get an estimate of those xi values:

xx       = normrnd(estims_xi(1,:),se_xi(1,:),J-1,1);                                            
mm       = normrnd(estims_xi(2,:),se_xi(2,:),M,1);                           
tt       = normrnd(estims_xi(3,:),se_xi(3,:),T,1);                                            
noise    = normrnd(0,0.5,J-1,M,T);                                     

xi_jmt_estimreconstructed  = zeros(J-1,M,T);

for   j = 1:J-1
  for market = 1:M
    for time = 1:T
        xi_jmt_estimreconstructed(j,market,time) = xx(j,1) + mm(market,1) + tt(time,1) + noise(j,market,time);
    end
  end
end

xi_jmt_estimreconstructed               = xi_jmt_estimreconstructed(:);                                        
xi_jmt_estimreconstructed               = xi_jmt_estimreconstructed(2:end,:);

first_xi_regression  = xi_jmt_estimreconstructed(first_index,:);
second_xi_regression = xi_jmt_estimreconstructed(second_index,:);

first_xi_R  = xi_jmt(first_index,:);
second_xi_R = xi_jmt(second_index,:);

first_xi_M  = mean(first_xi);
second_xi_M = mean(second_xi);

beta_est  =  X_J(1:4,:);
alpha_est = -X_J(5,:);
sigma_est =  X_J(6,:);
gamma_est =  X_J(7:end,:);
mc_est    = [[1 w_111]*gamma_est+base_omega_jmt(1,:); mc_J(1:J-2)];  




%% Herein, I start solving for new eq. prices & shares for first and second xi specifications:

inits_xi = ones(J-1,1);

for k = 1: size(first_xi,1)
    
solveforprices_first_xi    = @(pppp) (pppp -shares_RC_counterfactual(pppp,first_xi(k,:)) - mc_est);
    
tic
p_jmt_first_xi(k,:)                = fsolve(solveforprices_first_xi,inits_xi,opt);             
toc

s_jmt_first_xi(k,:) = sharescalculator_RC_counterfactual(p_jmt_first_xi(k,:)',first_xi(k,:));

end

for k = 1: size(second_xi,1)
    
solveforprices_second_xi    = @(pppp) (pppp -shares_RC_counterfactual(pppp,second_xi(k,:)) - mc_est);
    
tic
p_jmt_second_xi(k,:)                = fsolve(solveforprices_second_xi,inits_xi,opt);             
toc

s_jmt_second_xi(k,:) = sharescalculator_RC_counterfactual(p_jmt_second_xi(k,:)',second_xi(k,:));

end

%% Now calculate the new equilibrium prices and shares for the MEAN of those xi specifications:

solveforprices_first_xi_M    = @(pppp) (pppp -shares_RC_counterfactual(pppp,first_xi_M) - mc_est);
    
tic
p_jmt_first_xi_M          = fsolve(solveforprices_first_xi_M,inits_xi,opt);             
toc

s_jmt_first_xi_M = sharescalculator_RC_counterfactual(p_jmt_first_xi_M,first_xi_M);


solveforprices_second_xi_M    = @(pppp) (pppp -shares_RC_counterfactual(pppp,second_xi_M) - mc_est);
    
tic
p_jmt_second_xi_M          = fsolve(solveforprices_second_xi_M,inits_xi,opt);             
toc

s_jmt_second_xi_M = sharescalculator_RC_counterfactual(p_jmt_first_xi_M,second_xi_M);

%% Now I use the regression values as an estimate of first xi and second xi:

for k = 1: size(first_xi_regression,1)
    
solveforprices_first_xi_regression    = @(pppp) (pppp -shares_RC_counterfactual(pppp,first_xi_regression(k,:)) - mc_est);
    
tic
p_jmt_first_xi_regression(k,:)                = fsolve(solveforprices_first_xi_regression,inits_xi,opt);             
toc

s_jmt_first_xi_regression(k,:) = sharescalculator_RC_counterfactual(p_jmt_first_xi_regression(k,:)',first_xi_regression(k,:));

end

for k = 1: size(second_xi_regression,1)
    
solveforprices_second_xi_regression    = @(pppp) (pppp -shares_RC_counterfactual(pppp,second_xi_regression(k,:)) - mc_est);
    
tic
p_jmt_second_xi_regression(k,:)        = fsolve(solveforprices_second_xi_regression,inits_xi,opt);             
toc

s_jmt_second_xi_regression(k,:) = sharescalculator_RC_counterfactual(p_jmt_second_xi_regression(k,:)',second_xi_regression(k,:));

end

%% Lastly, calculate those profits & shares for the REAL values of xi:

for k = 1: size(first_xi_R,1)
    
solveforprices_first_xi_R    = @(pppp) (pppp -shares_RC_counterfactual(pppp,first_xi_R(k,:)) - mc_est);
    
tic
p_jmt_first_xi_R(k,:)                = fsolve(solveforprices_first_xi_R,inits_xi,opt);             
toc

s_jmt_first_xi_R(k,:) = sharescalculator_RC_counterfactual(p_jmt_first_xi_R(k,:)',first_xi_R(k,:));

end

for k = 1: size(second_xi_R,1)
    
solveforprices_second_xi_R    = @(pppp) (pppp -shares_RC_counterfactual(pppp,second_xi_R(k,:)) - mc_est);
    
tic
p_jmt_second_xi_R(k,:)                = fsolve(solveforprices_second_xi_R,inits_xi,opt);             
toc

s_jmt_second_xi_R(k,:) = sharescalculator_RC_counterfactual(p_jmt_second_xi_R(k,:)',second_xi_R(k,:));

end

%% Calculation of profits for new entrant firm - rivals: IN TOTAL THERE WILL BE 8 DIFFERENT ONES:

firstbrand_profits_baseline  = 0;                                          % The good is missing at baseline so it's 0
rivals_profits_baseline      = sum((p_jmt(1:J-2,:) - mc_J(1:J-2,:)).*s_jmt(1:J-2,:));  
                                                                           % This gives the sum of rival firm profits for M=1&T=1 in baseline specification.

                                                
mc_estim                     = mc_est';                                    % Transpozing for dimensional matching

profits_first_xi             = (p_jmt_first_xi - mc_estim).*s_jmt_first_xi;% This gives the profits for all 5 goods in first counterfactual
firstbrand_profits_first_xi  = mean(profits_first_xi(:,1));                % This gives the MEAN of the first brand for M=1&T=1 in counterfactual
                                                                           % specification for the first xi
rivals_profits_first_xi      = mean(sum(profits_first_xi(:,2:end),2));     % This gives the MEAN of the sum of rival firm profits for M=1&T=1 in 
                                                                           % counterfactual specification for the first xi

profits_second_xi            = (p_jmt_second_xi - mc_estim).*s_jmt_second_xi;     
                                                                           % This gives the profits for all 5 goods in second counterfactual

firstbrand_profits_second_xi = mean(profits_second_xi(:,1));               % This gives the MEAN of the first brand for M=1&T=1 in counterfactual
                                                                           % specification for the second xi
rivals_profits_second_xi     = mean(sum(profits_second_xi(:,2:end),2));    % This gives the MEAN of the sum of rival firm profits for M=1&T=1 in 
                                                                           % counterfactual specification for the second xi

                                                                           
profits_first_xi_M            = (p_jmt_first_xi_M - mc_est).*s_jmt_first_xi_M;
firstbrand_profits_first_xi_M = profits_first_xi_M(1,:);
rivals_profits_first_xi_M     = sum(profits_first_xi_M(2:end,:));


profits_second_xi_M            = (p_jmt_second_xi_M - mc_est).*s_jmt_second_xi_M;
firstbrand_profits_second_xi_M = profits_second_xi_M(1,:);
rivals_profits_second_xi_M     = sum(profits_second_xi_M(2:end,:));


profits_first_xi_regression             = (p_jmt_first_xi_regression - mc_estim).*s_jmt_first_xi_regression;% This gives the profits for all 5 goods in first counterfactual
firstbrand_profits_first_xi_regression  = mean(profits_first_xi_regression(:,1));   % This gives the MEAN of the first brand for M=1&T=1 in counterfactual
                                                                           % specification for the first xi "REGRESSION"
rivals_profits_first_xi_regression      = mean(sum(profits_first_xi_regression(:,2:end),2));     % This gives the MEAN of the sum of rival firm profits for M=1&T=1 in 
                                                                           % counterfactual specification for the first xi "REGRESSION"

profits_second_xi_regression            = (p_jmt_second_xi_regression - mc_estim).*s_jmt_second_xi_regression;     
                                                                           % This gives the profits for all 5 goods in second counterfactual

firstbrand_profits_second_xi_regression = mean(profits_second_xi_regression(:,1));               % This gives the MEAN of the first brand for M=1&T=1 in counterfactual
                                                                           % specification for the second xi "REGRESSION"
rivals_profits_second_xi_regression     = mean(sum(profits_second_xi_regression(:,2:end),2));    % This gives the MEAN of the sum of rival firm profits for M=1&T=1 in 
                                                                           % counterfactual specification for the second xi "REGRESSION"


profits_first_xi_R             = (p_jmt_first_xi_R - mc_estim).*s_jmt_first_xi_R;      
                                                                           % This gives the profits for all 5 goods in first counterfactual
                                                                           
firstbrand_profits_first_xi_R  = mean(profits_first_xi_R(:,1));            % This gives the MEAN of the first brand for M=1&T=1 in counterfactual
                                                                           % specification for the first xi (REAL ONE)
                                                                           
rivals_profits_first_xi_R      = mean(sum(profits_first_xi_R(:,2:end),2)); % This gives the MEAN of the sum of rival firm profits for M=1&T=1 in 
                                                                           % counterfactual specification for the first xi (REAL ONE)

profits_second_xi_R            = (p_jmt_second_xi_R - mc_estim).*s_jmt_second_xi_R;   
                                                                           % This gives the profits for all 5 goods in second counterfactual

firstbrand_profits_second_xi_R = mean(profits_second_xi_R(:,1));           % This gives the MEAN of the first brand for M=1&T=1 in counterfactual
                                                                           % specification for the second xi (REAL ONE)
                                                                           
rivals_profits_second_xi_R     = mean(sum(profits_second_xi_R(:,2:end),2));% This gives the MEAN of the sum of rival firm profits for M=1&T=1 in 
                                                                           % counterfactual specification for the second xi (REAL ONE)

%% Calculation for consumer surplus:

CS_baseline = zeros(J-2,I);
for i = 1:I
        CS_baseline(1:J-2,i) = [ones(J-2,1) x_jmt(1:J-2,:)] * beta_est ...
    - alpha_est .* p_jmt(1:J-2,:) + xi_J(1:J-2,:) - p_jmt(1:J-2,:)*sigma_est*heterogeneity(1,i);                       

end
CS_baseline = exp(CS_baseline);
CS_baseline = sum(CS_baseline,1);
CS_baseline = log(CS_baseline);
CS_baseline = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_baseline;
CS_baseline = mean(CS_baseline);   % Negative -> C is negative

CS_first_xi = zeros(J-1,I,size(first_xi,1));
for i = 1:I
 for u = 1:size(first_xi,1)
CS_first_xi(1:J-1,i,u) =  [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_first_xi(u,1:J-1)' + [first_xi(u,1);xi_J(1:J-2,:)] - p_jmt_first_xi(u,1:J-1)'*sigma_est*heterogeneity(1,i);      
 end
end
CS_first_xi = exp(CS_first_xi);
CS_first_xi = sum(CS_first_xi,1);
CS_first_xi = log(CS_first_xi);
CS_first_xi = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_first_xi;
CS_first_xi = mean(CS_first_xi);
CS_first_xi = mean(CS_first_xi);
CS_first_xi = CS_first_xi - CS_baseline;

CS_second_xi = zeros(J-1,I,size(second_xi,1));
for i = 1:I
 for u = 1:size(second_xi,1)
CS_second_xi(1:J-1,i,u) =  [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_second_xi(u,1:J-1)' + [second_xi(u,1);xi_J(1:J-2,:)] - p_jmt_second_xi(u,1:J-1)'*sigma_est*heterogeneity(1,i);      
 end
end
CS_second_xi = exp(CS_second_xi);
CS_second_xi = sum(CS_second_xi,1);
CS_second_xi = log(CS_second_xi);
CS_second_xi = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_second_xi;
CS_second_xi = mean(CS_second_xi);
CS_second_xi = mean(CS_second_xi);
CS_second_xi = CS_second_xi - CS_baseline;

CS_first_xi_M = zeros(J-1,I);
for i = 1:I
        CS_first_xi_M(1:J-1,i) = [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_first_xi_M(1:J-1,:) + [first_xi_M;xi_J(1:J-2,:)] - p_jmt_first_xi_M(1:J-1,:)*sigma_est*heterogeneity(1,i);                       
end
CS_first_xi_M = exp(CS_first_xi_M);
CS_first_xi_M = sum(CS_first_xi_M,1);
CS_first_xi_M = log(CS_first_xi_M);
CS_first_xi_M = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_first_xi_M;
CS_first_xi_M = mean(CS_first_xi_M); 
CS_first_xi_M = CS_first_xi_M - CS_baseline;

CS_second_xi_M = zeros(J-1,I);
for i = 1:I
        CS_second_xi_M(1:J-1,i) = [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_second_xi_M(1:J-1,:) + [second_xi_M;xi_J(1:J-2,:)] - p_jmt_second_xi_M(1:J-1,:)*sigma_est*heterogeneity(1,i);                       
end
CS_second_xi_M = exp(CS_second_xi_M);
CS_second_xi_M = sum(CS_second_xi_M,1);
CS_second_xi_M = log(CS_second_xi_M);
CS_second_xi_M = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_second_xi_M;
CS_second_xi_M = mean(CS_second_xi_M); 
CS_second_xi_M = CS_second_xi_M - CS_baseline;

CS_first_xi_regression = zeros(J-1,I,size(first_xi_regression,1));
for i = 1:I
 for u = 1:size(first_xi_regression,1)
CS_first_xi_regression(1:J-1,i,u) =  [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_first_xi_regression(u,1:J-1)' + [first_xi_regression(u,1);xi_J(1:J-2,:)] - p_jmt_first_xi_regression(u,1:J-1)'*sigma_est*heterogeneity(1,i);      
 end
end
CS_first_xi_regression = exp(CS_first_xi_regression);
CS_first_xi_regression = sum(CS_first_xi_regression,1);
CS_first_xi_regression = log(CS_first_xi_regression);
CS_first_xi_regression = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_first_xi_regression;
CS_first_xi_regression = mean(CS_first_xi_regression);
CS_first_xi_regression = mean(CS_first_xi_regression);
CS_first_xi_regression = CS_first_xi_regression - CS_baseline;

CS_second_xi_regression = zeros(J-1,I,size(second_xi_regression,1));
for i = 1:I
 for u = 1:size(second_xi_regression,1)
CS_second_xi_regression(1:J-1,i,u) =  [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_second_xi_regression(u,1:J-1)' + [second_xi_regression(u,1);xi_J(1:J-2,:)] - p_jmt_second_xi_regression(u,1:J-1)'*sigma_est*heterogeneity(1,i);      
 end
end
CS_second_xi_regression = exp(CS_second_xi_regression);
CS_second_xi_regression = sum(CS_second_xi_regression,1);
CS_second_xi_regression = log(CS_second_xi_regression);
CS_second_xi_regression = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_second_xi_regression;
CS_second_xi_regression = mean(CS_second_xi_regression);
CS_second_xi_regression = mean(CS_second_xi_regression);
CS_second_xi_regression = CS_second_xi_regression - CS_baseline;

CS_first_xi_R = zeros(J-1,I,size(first_xi_R,1));
for i = 1:I
 for u = 1:size(first_xi_R,1)
CS_first_xi_R(1:J-1,i,u) =  [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_first_xi_R(u,1:J-1)' + [first_xi_R(u,1);xi_J(1:J-2,:)] - p_jmt_first_xi_R(u,1:J-1)'*sigma_est*heterogeneity(1,i);      
 end
end
CS_first_xi_R = exp(CS_first_xi_R);
CS_first_xi_R = sum(CS_first_xi_R,1);
CS_first_xi_R = log(CS_first_xi_R);
CS_first_xi_R = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_first_xi_R;
CS_first_xi_R = mean(CS_first_xi_R);
CS_first_xi_R = mean(CS_first_xi_R);
CS_first_xi_R = CS_first_xi_R - CS_baseline;

CS_second_xi_R = zeros(J-1,I,size(second_xi_R,1));
for i = 1:I
 for u = 1:size(second_xi_R,1)
CS_second_xi_R(1:J-1,i,u) =  [ones(J-1,1) base_x_jmt(1:J-1,:)] * beta_est ...
    - alpha_est .* p_jmt_second_xi_R(u,1:J-1)' + [second_xi_R(u,1);xi_J(1:J-2,:)] - p_jmt_second_xi_R(u,1:J-1)'*sigma_est*heterogeneity(1,i);      
 end
end
CS_second_xi_R = exp(CS_second_xi_R);
CS_second_xi_R = sum(CS_second_xi_R,1);
CS_second_xi_R = log(CS_second_xi_R);
CS_second_xi_R = (1./(alpha_est + sigma_est.* heterogeneity)).*CS_second_xi_R;
CS_second_xi_R = mean(CS_second_xi_R);
CS_second_xi_R = mean(CS_second_xi_R);
CS_second_xi_R = CS_second_xi_R - CS_baseline;
