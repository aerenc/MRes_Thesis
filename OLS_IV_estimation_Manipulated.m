function  [m_OLSestimates,m_xi_jmt_error,m_Twoslsestimatessecondstage] = OLS_IV_estimation_Manipulated(Total,M,T,adj,m_index3,m_s_jmt,m_s0_mt,m_x_jmt,m_p_jmt,...
                                                                                                       m_w_jmt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we need to construct log(S_{jm}/S_{0m}) dependent variable:

m_dependentlogratio = zeros(Total-M*T-adj*M*T,1);

for i = 1 : M*T
    m_dependentlogratio(m_index3(i,1):m_index3(i,2),:) = log(m_s_jmt(m_index3(i,1):m_index3(i,2),:))-log(m_s0_mt(i,1));
end

m_X            = [ones(Total-M*T-adj*M*T,1) m_x_jmt m_p_jmt];   % stacking regressors together
m_OLSestimates = (m_X'*m_X)\m_X'*m_dependentlogratio;           % classical OLS formula giving us OLS estimates on aggregate Logit model

m_xi_jmt_error = m_dependentlogratio - m_X * m_OLSestimates ;   % This gives estimated unobserved product characteristics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF FIRST STAGE IV FOR MANIPULATED DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, to be able to remove the endogeneity bias in our model, we regress
% price additionally on cost shifters and collect the corresponding estimates for prices

m_Z = [ones(Total-M*T-adj*M*T,1) m_x_jmt m_w_jmt];                                 
m_price_hat = m_Z/(m_Z'*m_Z)*m_Z'*m_p_jmt;               % part of price that has no correlation with error term
m_Twoslsestimatesfirststage = (m_Z'*m_Z)\m_Z'*m_p_jmt;   % first stage estimates of coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL FOR MANIPULATED DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since we have the estimates for prices, we can now estimate our model by
% also estimating the second stage of IV regression:

m_X_hat = [ones(Total-M*T-adj*M*T,1) m_x_jmt m_price_hat];
m_Twoslsestimatessecondstage = (m_X_hat' * m_X_hat)\m_X_hat'*m_dependentlogratio;
