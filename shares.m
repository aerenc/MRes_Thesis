function [inside_s_jmt] = shares(p,x_jmt,xi_jmt,index,index2,index3,adj)

global beta alpha Total J M T

delta_jmt          = zeros(Total-adj*M*T,1);
delta_jmt(index,:) = 0;

for z = 1: M*T                                                                                             % prices
delta_jmt(index2(z,1):index2(z,2),:)  = [ones(J-adj-1,1) x_jmt(index3(z,1):index3(z,2),:)] * beta - alpha .* p(index3(z,1):index3(z,2),:) + xi_jmt(index3(z,1):index3(z,2),:);
end

nom               = exp(delta_jmt);
cumsum_delta_jmt  = cumsum(nom);
cumsum_1          = cumsum_delta_jmt(J-adj,:);
totalsum          = diff(cumsum_delta_jmt(index,:)); 
totalsum          = [cumsum_1;totalsum];         % M*T total share of each j \in J_{mt}
Totalsum          = repelem(totalsum,J-adj);     % stretching it out to (J*  M*T)*1 (since we have J-1 products in each market)
denom             = 1 + Totalsum;
s__jmt            = nom./denom;                  % getting ALL shares, i.e. including s0_mt 

inside_s_jmt      = zeros(Total-M*T-adj*M*T,1);  % now getting shares of inside goods

for i = 1: M*T
inside_s_jmt(index3(i,1):index3(i,2))  = s__jmt(index2(i,1):index2(i,2),1); % getting shares of inside goods
end

end