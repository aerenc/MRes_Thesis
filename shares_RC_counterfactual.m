function [inside_s_jmt] = shares_RC_counterfactual(p,base_x_jmt,xi_J,xi_vals)

global beta_est alpha_est J I heterogeneity sigma_est  

delta_jmt = zeros(J,I);

delta_jmt(J,:) = 0;

for i = 1:I
delta_jmt(1,i) =  [1 base_x_jmt(1,:)] * beta_est + alpha_est .* p(1,:) + xi_vals - p(1,:)*sigma_est*heterogeneity(1,i);
end

for i = 1:I
delta_jmt(2:J-1,i)  = [ones(J-2,1) base_x_jmt(2:J-1,:)] * beta_est ...
    + alpha_est .* p(2:end,:) + xi_J(1:J-2,:)...
  - p(2:end,:)*sigma_est*heterogeneity(1,i); 
end


nom               = exp(delta_jmt);
cums = cumsum(nom,1);
cums = cums(end,:);
denom             = cums;
s__jmt            = nom./denom;                                            

s___jmt           = mean(s__jmt,2);

inside_s_jmt      = zeros(J-1,1);                  % Now, getting shares of inside goods

inside_s_jmt(1:J-1,:)  = s___jmt(1:J-1,:);         % Getting shares of inside goods



end







