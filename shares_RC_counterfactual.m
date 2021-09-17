function [inside_s_jmt] = shares_RC_counterfactual(p,xi_vals)

global beta_est alpha_est J sigma_est base_x_jmt xi_J v ns

delta_mu_jmt = zeros(J,ns);

delta_mu_jmt(J,:) = 0;

for i = 1:ns
delta_mu_jmt(1,i) =  [1 base_x_jmt(1,:)] * beta_est - alpha_est .* p(1,:) + xi_vals - p(1,:)*sigma_est*v(1,i);
end

for i = 1:ns
delta_mu_jmt(2:J-1,i)  = [ones(J-2,1) base_x_jmt(2:J-1,:)] * beta_est ...
    - alpha_est .* p(2:end,:) + xi_J(1:J-2,:)...
  - p(2:end,:)*sigma_est*v(1,i); 
end


nom               = exp(delta_mu_jmt);
cums = cumsum(nom,1);
cums = cums(end,:);
denom             = cums;
s_ijmt            = nom./denom;                                            

nom_general       = mean(s_ijmt,2);
denom_general     = mean((alpha_est+sigma_est.*v).*s_ijmt.*(1-s_ijmt),2);
share_general     = nom_general./denom_general; 

inside_s_jmt      = zeros(J-1,1);                  

inside_s_jmt(1:J-1,:)  = share_general(1:J-1,:);         % Getting shares of inside goods



end







