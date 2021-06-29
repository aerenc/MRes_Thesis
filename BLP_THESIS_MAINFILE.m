% BLP THESIS
% SEVERAL APPROACHES TO UNOBSERVED PRODUCT CHARACTERISTICS WITH SIMULATED DATA UNDER LOGIT & RC SPECIFICATION
% Ali Eren Camur
% Tilburg University - CentER

clc;
clear;

global x_jmt beta alpha xi_jmt Total J M T index index2 index3 inits adj m_x_jmt m_xi_jmt m_index m_index2 m_index3  

rng('default');

%% SIMULATION SETTING - LOGIT

% Generating the data for logit case first:

alpha   = 2;                    % true alpha
beta    = [3,3,0.5,0.5]';       % true betas
gamma   = [5,1,0.5,0.5]';       % true gamma
theta   = [alpha; beta; gamma]; % true theta

I       = 200;                  % Total number of consumers in each market
J       = 20;                   % Initial total good + outside share in each market
M       = 20;                   % Total markets
T       = 10;                   % Total time period
Total   = J*M*T;                % Total goods present a priori, INCLUDING outside good

%% INDEXES:

index   = [J:J:Total]';             % indexing each good in every market M*T combination
indexx  = [J-1:J-1:Total-M*T]';     % indexing that will be used in calculating outside good share in this setting

index2             = zeros(M*T,2);  % indexing the inside goods in each market in matrix form
index2(1,1)        = 1;
index2(1,2)        = J-1;
for i = 2 : M*T
    index2(i,1)    = index(i-1,1)+1;       
    index2(i,2)    = index(i,1)-1; 
end
index3             = zeros(M*T,2);  % indexing optimizer prices and shares etc.
index3(1,1)        = 1;
index3(1,2)        = J-1;
for i = 2 : M*T
   index3(i,1)     = index3(i-1,2)+1;
   index3(i,2)     = index3(i,1)+J-2;
end

adj       = 5;                                        % 5 goods will be removed from each market & time period randomly as in manipulated case

m_index   = [J-adj:J-adj:Total-adj*M*T]';             % indexing each good in every market M*T combination for manipulated case

m_indexx  = [J-1-adj:J-1-adj:Total-M*T-adj*M*T]';     % indexing that will be used in calculating outside good share for manipulated case

m_index2             = zeros(M*T,2);                  % indexing the inside goods in each market in matrix form for manipulated case
m_index2(1,1)        = 1;
m_index2(1,2)        = J-1-adj;

for i = 2 : M*T
    m_index2(i,1)    = m_index(i-1,1)+1;       
    m_index2(i,2)    = m_index(i,1)-1; 
end

m_index3             = zeros(M*T,2);                  % indexing optimizer prices for manipulated case
m_index3(1,1)        = 1;
m_index3(1,2)        = J-adj-1;

for i = 2 : M*T
   m_index3(i,1)     = m_index3(i-1,2)+1;
   m_index3(i,2)     = m_index3(i,1) + J-adj-2;
end

mm_drop_index = zeros(adj,M);     % these elements will be dropped in manipulation: the dropped brands will be same for all time periods for each market

for i = 1 : M
mm_drop_index(:,i) = randperm(J-1,adj)';
end

mm_drop_index = [16,18,6,13,4,14,15,16,11,12,5,13,12,10,9,9,7,2,12,15;12,16,12,16,13,8,1,11,10,3,18,5,2,11,2,7,9,5,17,7;2,6,13,3,14,9,10,10,1,5,1,16,15,19,3,5,18,1,16,6;6,17,18,1,16,6,14,8,18,15,8,15,7,8,15,14,8,16,11,13;19,5,10,17,12,17,2,1,8,14,13,12,19,1,11,6,5,3,18,11];
mm_drop_index = sort(mm_drop_index);

m_drop_index  = zeros(adj,M*T);

t_index        = zeros(M,2);
t_index(1,1)   = 1;
t_index(1,2)   = T;

for  i = 2:M
   t_index(i,1) = t_index(i-1,2) + 1;
   t_index(i,2) = t_index(i,1) + T - 1;
end

for i = 1 : M
    m_drop_index(:,i:M:M*T-M+i) = repmat(mm_drop_index(:,i),1,T);
end

%% SIMULATING THE DATA

% excluding the outside good <=> Total - M*T

x_jmt_1 = unifrnd(0,1,Total-M*T,1);
x_jmt_2 = unifrnd(0,1,Total-M*T,1);
x_jmt_3 = unifrnd(0,1,Total-M*T,1);
x_jmt   = [x_jmt_1 x_jmt_2 x_jmt_3];  % Generated observed product characteristics

w_jmt_1 = unifrnd(0,1,Total-M*T,1);
w_jmt_2 = unifrnd(0,1,Total-M*T,1);
w_jmt_3 = unifrnd(0,1,Total-M*T,1);
w_jmt   = [w_jmt_1 w_jmt_2 w_jmt_3];  % Generated cost shifters

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

xi_jmt               = xi_jmt(:);                                     % unobserved product characteristics vector
omega_jmt            = normrnd(0,1,Total-M*T,1);                      % unobserved costs
mc_jmt               = [ones(Total-M*T,1) w_jmt] * gamma + omega_jmt; % marginal cost specification
inits                = ones(Total-M*T,1);
solveforprices       = @(ppp) (ppp -(1./(alpha.*(1 - s_jmt(ppp)))) - mc_jmt);
%opts                = optimset('Algorithm ','interior-point ','Display ','iter ');
opts                 = optimset('Display','off');
%opts                = optimset('Display','iter','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-10,'MaxFunEvals',1000000000);

tic
p_jmt                = fsolve(solveforprices,inits,opts); % this gives us prices
toc

s_jmt                = s_jmt(p_jmt);                      % gives us inner market shares
cumsum_s_jmt         = cumsum(s_jmt);
cumsum_s_jmt_1       = cumsum_s_jmt(J-1,:);
marketwisesum_s_jmt  = diff(cumsum_s_jmt(indexx,:)); 
marketwisesum_s_jmt  = [cumsum_s_jmt_1;marketwisesum_s_jmt]; 
s0_mt                = 1 - marketwisesum_s_jmt;           % getting outside shares for each M*T combo

%% UNMANIPULATED DATA IN 3D-ARRAY ((J-1)*M*T) FORM:

array_s_jmt     = reshape(s_jmt,J-1,M,T);         % market shares
array_p_jmt     = reshape(p_jmt,J-1,M,T);         % prices 
array_xi_jmt    = reshape(xi_jmt,J-1,M,T);        % unobs. prod. characteristics
array_omega_jmt = reshape(omega_jmt,J-1,M,T);     % unobs. costs
array_mc_jmt    = reshape(mc_jmt,J-1,M,T);        % marginal costs

array_x_jmt     = zeros(J-1,M,T,size(x_jmt,2));   % observed prod. characteristics
x11_jmt = x_jmt(:,1);
x12_jmt = x_jmt(:,2);
x13_jmt = x_jmt(:,3);
array_x_jmt(:,:,:,1) = reshape(x11_jmt,J-1,M,T);
array_x_jmt(:,:,:,2) = reshape(x12_jmt,J-1,M,T);  
array_x_jmt(:,:,:,3) = reshape(x13_jmt,J-1,M,T);  
array_w_jmt = zeros(J-1,M,T,size(w_jmt,2));       % observed cost shifters
w11_jmt = w_jmt(:,1);
w12_jmt = w_jmt(:,2);
w13_jmt = w_jmt(:,3);
array_w_jmt(:,:,:,1) = reshape(w11_jmt,J-1,M,T);  
array_w_jmt(:,:,:,2) = reshape(w12_jmt,J-1,M,T);  
array_w_jmt(:,:,:,3) = reshape(w13_jmt,J-1,M,T);  

%% CONSTRUCTING THE BLP INSTRUMENTS FOR UNMANIPULATED DATA: 

% % Note that we will have two types of instruments: own product characteristics and characteristics sum of rivals' characteristics across markets and time periods
% 
% cumsum_z1_jmt     = cumsum(x_jmt);
% cumsum_z1         = cumsum_z1_jmt(J-1,:);
% totalsum_z1       = diff(cumsum_z1_jmt(indexx,:)); 
% totalsum_z1       = [cumsum_z1;totalsum_z1]; %(M*T)*4 -> contains the sum of product characteristics for every M*T combination
% 
% for  i = 1: M*T
%      rival_jmt(index3(i,1):index3(i,2),:) = totalsum_z1(i,:) - x_jmt(index3(i,1):index3(i,2),:);   % sum of rivals' observed product characteristics 
% end
% 
% Z_blp = [ones(Total-M*T,1) x_jmt rival_jmt x_jmt.^2 rival_jmt.^2]; %BLP instruments -> contains own & rivals prod. characteristics and their squares.


% %% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL 

% % Now we need to construct log(S_{jm}/S_{0m}) dependent variable:
% 
% dependentlogratio = zeros(Total-M*T,1);
% 
% for i = 1 : M*T
%     dependentlogratio(index3(i,1):index3(i,2),:) = log(s_jmt(index3(i,1):index3(i,2),:))-log(s0_mt(i,1));
% end
% 
% X            = [ones(Total-M*T,1) x_jmt p_jmt];              % stacking regressors together
% OLSestimates = (X'*X)\X'*dependentlogratio; % classical OLS formula giving us OLS estimates on aggregate Logit model
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% OLS ESTIMATION OF FIRST STAGE IV %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % First, to be able to remove the endogeneity bias in our model, we regress
% % price additionally on cost shifters and collect the corresponding estimates for prices
% 
% Z = [ones(Total-M*T,1) w_jmt x_jmt];                                 
% price_hat = Z/(Z'*Z)*Z'*p_jmt;               % part of price that has no correlation with error term
% Twoslsestimatesfirststage = (Z'*Z)\Z'*p_jmt; % First stage estimates of coefficients
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Since we have the estimates for prices, we can now estimate our model by
% % also estimating the second stage of IV regression:
% 
% X_hat = [ones(Total-M*T,1) x_jmt price_hat];
% Twoslsestimatessecondstage = (X_hat' * X_hat)\X_hat'*dependentlogratio;

%% MANIPULATION BEGINS

% Bir onceki manipulation'da her market ve time period'dan random bir
% sekilde 5 tane eleman eledim. Ama given bir markette bu good t'de exist
% edip t+1'de etmeyip t+2'de ederse (ki onceki manipulation bunu mumkun
% kiliyor) bu intuitive olarak pek makul olmayabilir. Simdi her M icin tum
% time periodlarda ayni elemanlari dusurerek datayi bastan manipule
% edecegiz. Her M icin adj = 5 devam ediyor

% We reshape the data into (J-1) * (M*T) matrices for each element of
% x and w:

reshaped_x1_jmt    = reshape(x_jmt(:,1),J-1,M*T); 
reshaped_x2_jmt    = reshape(x_jmt(:,2),J-1,M*T); 
reshaped_x3_jmt    = reshape(x_jmt(:,3),J-1,M*T); 
reshaped_w1_jmt    = reshape(w_jmt(:,1),J-1,M*T); 
reshaped_w2_jmt    = reshape(w_jmt(:,2),J-1,M*T); 
reshaped_w3_jmt    = reshape(w_jmt(:,3),J-1,M*T); 

reshaped_xi_jmt    = reshape(xi_jmt,J-1,M*T); 
reshaped_omega_jmt = reshape(omega_jmt,J-1,M*T);

% Now dropping the elements from unmanipulated case:

for  i = 1 : M*T
    reshaped_x1_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_x2_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_x3_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_w1_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_w2_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_w3_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_xi_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    reshaped_omega_jmt(m_drop_index(:,i),i) = 0;
end

%%
% m stands for manipulated:

m_x1_jmt = reshaped_x1_jmt(:);
m_x1_jmt = m_x1_jmt(m_x1_jmt~=0);
m_x2_jmt = reshaped_x2_jmt(:);
m_x2_jmt = m_x2_jmt(m_x2_jmt~=0);
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

m_inits               = ones(Total-M*T-adj*M*T,1);
m_solveforprices      = @(pp) (pp -(1./(alpha.*(1 - m_s_jmt(pp)))) - m_mc_jmt);

tic
m_p_jmt               = fsolve(m_solveforprices,m_inits,opts); % this gives us prices
toc

m_s_jmt               = m_s_jmt(m_p_jmt);          % gives us inner market shares
cumsum_m_s_jmt        = cumsum(m_s_jmt);
cumsum_m_s_jmt_1      = cumsum_m_s_jmt(J-adj-1,:);
marketwisesum_m_s_jmt = diff(cumsum_m_s_jmt(m_indexx,:)); 
marketwisesum_m_s_jmt = [cumsum_m_s_jmt_1;marketwisesum_m_s_jmt]; 
m_s0_mt               = 1 - marketwisesum_m_s_jmt;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we need to construct log(S_{jm}/S_{0m}) dependent variable:

m_dependentlogratio = zeros(Total-M*T-adj*M*T,1);

for i = 1 : M*T
    m_dependentlogratio(m_index3(i,1):m_index3(i,2),:) = log(m_s_jmt(m_index3(i,1):m_index3(i,2),:))-log(m_s0_mt(i,1));
end

m_X            = [ones(Total-M*T-adj*M*T,1) m_x_jmt m_p_jmt];              % stacking regressors together
m_OLSestimates = (m_X'*m_X)\m_X'*m_dependentlogratio;       % classical OLS formula giving us OLS estimates on aggregate Logit model

m_xi_jmt_error = m_dependentlogratio - m_X * m_OLSestimates ;  % This gives estimated unobserved product characteristics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF FIRST STAGE IV %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, to be able to remove the endogeneity bias in our model, we regress
% price additionally on cost shifters and collect the corresponding estimates for prices

m_Z = [ones(Total-M*T-adj*M*T,1) m_x_jmt m_w_jmt];                                 
m_price_hat = m_Z/(m_Z'*m_Z)*m_Z'*m_p_jmt;               % part of price that has no correlation with error term
m_Twoslsestimatesfirststage = (m_Z'*m_Z)\m_Z'*m_p_jmt; % First stage estimates of coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since we have the estimates for prices, we can now estimate our model by
% also estimating the second stage of IV regression:

m_X_hat = [ones(Total-M*T-adj*M*T,1) m_x_jmt m_price_hat];
m_Twoslsestimatessecondstage = (m_X_hat' * m_X_hat)\m_X_hat'*m_dependentlogratio;

%%









