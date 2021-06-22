function [s__jmt] = s_jmt(p)

global x_jmt beta alpha xi_jmt Total J


delta_jmt         = x_jmt * beta - alpha .* p + xi_jmt; 
nom               = exp(delta_jmt);
cumsum_delta_jmt  = cumsum(nom);
cumsum_1          = cumsum_delta_jmt(J,:);
index             = [J:J:Total]';
totalsum          = diff(cumsum_delta_jmt(index,:)); 
totalsum          = [cumsum_1;totalsum]; % M*T total share of each j \in J_{mt}
Totalsum          = repelem(totalsum,J); % stretching it out to J*M*T)*1 (since we have 15 products in each market a priori)
denom             = 1 + Totalsum;
s__jmt            = nom./denom;














end



% p         = rand(Total,1);
% delta_mt  = x_jmt * beta + alpha .* p + xi_jmt;
% nom       = exp(delta_mt);
%  
% cumsum_delta_jmt  = cumsum(nom);
% 
% 
% 
% cumsum_1  = cumsum_delta_jmt(15,:);
% index     = [15:15:3000]';
% totalsum  = diff(cumsum_delta_jmt(index,:));
% totalsum  = [cumsum_1;totalsum];
% Totalsum  = repelem(totalsum,15); % \sum blabla
% denom     = ones(Total,1) + Totalsum;
% s_jmt     = nom./denom;
% 
% 
% solveforprices = @(p) (log(p - (ones(Total,1)./(alpha.*(ones(Total,1) - s_jmt)))) - w_jmt * gamma - omega_jmt);
% 
% %solveforprices = @(p) (log(p - (1./(alpha.*(1 - s_mt./(1+ repelem([cumsum(exp(x_jmt * beta + p .* alpha + xi_jmt))(15,:); diff(cumsum(exp(x_jmt*beta+p.*alpha+xi_jmt))([15:15:3000]',:))],15))))))) - w_jmt * gamma - omega_jmt);
% 
% 
% 
% %opts      = optimset('Display','off');
% opts    = optimset('Display','iter','TolCon',1E-14,'TolFun',1E-14,'TolX',1E-14);
% 
% tic
% prices  = fsolve(solveforprices,rand(Total,1),opts);
% toc
% 
% 
% 
% 

%solveforprices = @(p) (log(p - (1./(alpha.*(1 - exp(x_jmt * beta + alpha .* p + xi_jmt)./ ...
%    (1 + repelem([[[cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(15,:)];...
%    [diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:))]]],15)))))) - w_jmt * gamma - omega_jmt);


















%%

% delta_mt  = x_jmt * beta + alpha .* p + xi_jmt;
% nom       = exp(x_jmt * beta + alpha .* p + xi_jmt);
%  
% cumsum_delta_jmt  = cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt));
% cumsum_1  = cumsum_delta_jmt(15,:);
% 
% 
% 
% index     = [15:15:3000]';
% 
% 
% totalsum  = diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:));
% totalsum  = [cumsum_1;diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:))];
% Totalsum  = repelem([cumsum_1;diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:))],15); % \sum blabla
% denom     = 1 + repelem([cumsum_1;diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:))],15);
% s_jmt     = (exp(x_jmt * beta + alpha .* p + xi_jmt))./(1 + repelem([cumsum_1;diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:))],15));
% 
% 
% solveforprices = @(p) (log(p - (ones(Total,1)./(alpha.*(ones(Total,1) - (exp(x_jmt * beta + alpha .* p + xi_jmt))./(1 + repelem([cumsum_1;diff(cumsum(exp(x_jmt * beta + alpha .* p + xi_jmt))(index,:))],15)))))) - w_jmt * gamma - omega_jmt);
% 
% 
% 
% %opts      = optimset('Display','off');
% opts    = optimset('Display','iter','TolCon',1E-14,'TolFun',1E-14,'TolX',1E-14);
% 
% tic
% prices  = fsolve(solveforprices,rand(Total,1),opts);
% toc
% 




