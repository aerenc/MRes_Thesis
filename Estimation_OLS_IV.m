function  [OLSestimates,xi_jmt_error,Twoslsestimatessecondstage] = Estimation_OLS_IV(Total,M,T,one,index3,s_jmt,s0_mt,x_jmt,p_jmt,w_jmt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR "REAL" DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we need to construct log(S_{jm}/S_{0m}) dependent variable (i.e. y):

y            = zeros(Total-M*T-one,1);

for i = 1 : M*T
    y(index3(i,1):index3(i,2),:) = log(s_jmt(index3(i,1):index3(i,2),:))-log(s0_mt(i,1));
end

m_X          = [ones(Total-M*T-one,1) x_jmt p_jmt];                        % Stacking regressors together
OLSestimates = (m_X'*m_X)\m_X'*y;                                          % Classical OLS formula giving us OLS estimates on aggregate Logit model

xi_jmt_error = y - m_X * OLSestimates ;                                    % This gives estimated unobserved product characteristics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF FIRST STAGE IV FOR "REAL" DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, to be able to remove the endogeneity bias in our model, we regress
% price additionally on cost shifters and collect the corresponding estimates for prices

m_Z                         = [ones(Total-M*T-one,1) x_jmt w_jmt];                                 
m_price_hat                 = m_Z/(m_Z'*m_Z)*m_Z'*p_jmt;                   % Part of price that has no correlation with error term
m_Twoslsestimatesfirststage = (m_Z'*m_Z)\m_Z'*p_jmt;                       % First stage estimates of coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since we have the estimates for prices, we can now estimate our model by
% also estimating the second stage of IV regression:

m_X_hat                    = [ones(Total-M*T-one,1) x_jmt m_price_hat];
Twoslsestimatessecondstage = (m_X_hat' * m_X_hat)\m_X_hat'*y;
