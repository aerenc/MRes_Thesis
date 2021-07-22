function [inside_s_jmt] = shares_RC(p,x_jmt,xi_jmt,sigmaa,index,index2,index3)

global beta alpha Total J M T I heterogeneity
% p   = normrnd(5,1,Total-M*T-1,1);

delta_jmt          = zeros(Total-1,I);                                        % This will be filled with mean market values
delta_jmt(index,:) = 0;

for i = 1:I
delta_jmt(index2(1,1):index2(1,2),i)  = [ones(J-2,1) x_jmt(index3(1,1):index3(1,2),:)] * beta ...
                                                     - alpha .* p(index3(1,1):index3(1,2),:) + xi_jmt(index3(1,1):index3(1,2),:)...
 + repmat(sigmaa(1),J-2,1).*heterogeneity(1,i) - p(index3(1,1):index3(1,2),:)*sigmaa(2).*heterogeneity(2,i); 
end
                                                                     % Since this specific M*T
                                                                     % market lacks
                                                                     % one good
                                                                     % compared to
                                                                     % the other markets
                                                                      
for k = 2:M*T
    for i = 1:I
                delta_jmt(index2(k,1):index2(k,2),i)  = [ones(J-1,1) x_jmt(index3(k,1):index3(k,2),:)] * beta ...
                                                     - alpha .* p(index3(k,1):index3(k,2),:) + xi_jmt(index3(k,1):index3(k,2),:)...
  + repmat(sigmaa(1),J-1,1).*heterogeneity(1,i) - p(index3(k,1):index3(k,2),:).*sigmaa(2).*heterogeneity(2,i);   
    end
end



nom               = exp(delta_jmt);
cumsum_delta_jmt  = cumsum(nom,1);
cumsum_1          = cumsum_delta_jmt(J-1,:);
totalsum          = diff(cumsum_delta_jmt(index,:),1,1); 
totalsum          = [cumsum_1;totalsum];                                   % M*T total share of each j \in J_{mt}
Totalsum          = repelem(totalsum,J,1);                                 % Stretching it out to (J*M*T)*1 (since we have J-1 products in each market)
Totalsum          = Totalsum(2:end,:);
denom             = 1 + Totalsum;
s__jmt            = nom./denom;                                            % Getting ALL shares, i.e. including s0_mt 

s___jmt           = mean(s__jmt,2);

inside_s_jmt      = zeros(Total-M*T-1,1);                                  % Now, getting shares of inside goods

for i = 1: M*T
inside_s_jmt(index3(i,1):index3(i,2))  = s___jmt(index2(i,1):index2(i,2),1);% Getting shares of inside goods
end

end
