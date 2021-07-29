function [inside_shares] = sharescalculator_RC_counterfactual(p,xi_vals)

global beta_est alpha_est J I heterogeneity sigma_est base_x_jmt xi_J

delta_mu_jmt = zeros(J,I);

delta_mu_jmt(J,:) = 0;

for i = 1:I
delta_mu_jmt(1,i) =  [1 base_x_jmt(1,:)] * beta_est - alpha_est .* p(1,:) + xi_vals - p(1,:)*sigma_est*heterogeneity(1,i);
end

for i = 1:I
delta_mu_jmt(2:J-1,i)  = [ones(J-2,1) base_x_jmt(2:J-1,:)] * beta_est ...
    - alpha_est .* p(2:end,:) + xi_J(1:J-2,:)...
  - p(2:end,:)*sigma_est*heterogeneity(1,i); 
end


nom               = exp(delta_mu_jmt);
cums              = cumsum(nom,1);
cums              = cums(end,:);
denom             = cums;
s_ijmt            = nom./denom;                                            

s__jmt            = mean(s_ijmt,2);


inside_shares     = zeros(J-1,1);                   

inside_shares(1:J-1,:)  = s__jmt(1:J-1,:);          



end







