function [s___jmt] = m2_s_jmt(p)

global m2_x_jmt beta alpha m2_xi_jmt Total J M T m_index m_index2 m_index3 adj

delta_jmt          = zeros(Total-adj*M*T,1);
delta_jmt(m_index,:) = 0;
for z = 1: M*T
delta_jmt(m_index2(z,1):m_index2(z,2),:)  = [ones(J-adj-1,1) m2_x_jmt(m_index3(z,1):m_index3(z,2),:)] * beta - alpha .* p(m_index3(z,1):m_index3(z,2),:) + m2_xi_jmt(m_index3(z,1):m_index3(z,2),:);
end
nom               = exp(delta_jmt);
cumsum_delta_jmt  = cumsum(nom);
cumsum_1          = cumsum_delta_jmt(J-adj,:);
totalsum          = diff(cumsum_delta_jmt(m_index,:)); 
totalsum          = [cumsum_1;totalsum]; % M*T total share of each j \in J_{mt}
Totalsum          = repelem(totalsum,J-adj); % stretching it out to ((J-adj)*  M*T)*1 (since we have J-adj products in each market)
denom             = 1 + Totalsum;
s__jmt            = nom./denom;

s___jmt  = zeros(Total-M*T-adj*M*T,1);

for i = 1: M*T
s___jmt(m_index3(i,1):m_index3(i,2))  = s__jmt(m_index2(i,1):m_index2(i,2),1); % getting shares of inside goods
end

