function  [OLSestimates,xi_jmt_error_logit,Twoslsestimatessecondstage,xi_jmt_error_twosls,se_ols,se_twosls] = Estimation_OLS_IV(Total,M,T,one,index3,s_jmt,s0_mt,x_jmt,p_jmt,w_jmt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR "REAL" DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we need to construct log(S_{jm}/S_{0m}) dependent variable (i.e. y):

y            = zeros(Total-M*T-one,1);

for i = 1 : M*T
    y(index3(i,1):index3(i,2),:) = log(s_jmt(index3(i,1):index3(i,2),:))-log(s0_mt(i,1));
end

X            = [ones(Total-M*T-one,1) x_jmt p_jmt];                        % Stacking regressors together
OLSestimates = (X'*X)\X'*y;                                                % Classical OLS formula giving us OLS estimates on aggregate Logit model
xi_jmt_error_logit = y - X * OLSestimates ;                                      % This gives estimated unobserved product characteristics

n_logit=size(X,1);
k_logit=size(X,2);

covbeta_ols=(xi_jmt_error_logit'*xi_jmt_error_logit/(n_logit-k_logit))*inv(X'*X);
se_ols=sqrt(diag(covbeta_ols));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF FIRST STAGE IV FOR "REAL" DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, to be able to remove the endogeneity bias in our model, we regress
% price additionally on cost shifters and collect the corresponding estimates for prices

Z                          = [ones(Total-M*T-one,1) x_jmt w_jmt];                                 
price_hat                  = Z/(Z'*Z)*Z'*p_jmt;                            % Part of price that has no correlation with error term
Twoslsestimatesfirststage  = (Z'*Z)\Z'*p_jmt;                              % First stage estimates of coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since we have the estimates for prices, we can now estimate our model by
% also estimating the second stage of IV regression:

X_hat                      = [ones(Total-M*T-one,1) x_jmt price_hat];
Twoslsestimatessecondstage = (X_hat' * X_hat)\X_hat'*y;
xi_jmt_error_twosls = y - X_hat * Twoslsestimatessecondstage ;                                      % This gives estimated unobserved product characteristics

n_twosls=size(X,1);
k_twosls=size(X,2);

covbeta_twosls=(xi_jmt_error_twosls'*xi_jmt_error_twosls/(n_twosls-k_twosls))*inv(X_hat'*X_hat);
se_twosls=sqrt(diag(covbeta_twosls));