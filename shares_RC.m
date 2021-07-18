function [inside_s_jmt] = shares_RC(p,x_jmt,xi_jmt,sigma,index,index2,index3)

global beta alpha Total J M T

RC_exp          = zeros(Total-1,1);                                        % This will be filled with mean market values
RC_exp(index,:) = 0;

RC_exp(index2(1,1):index2(1,2),:)  = [ones(J-2,1) x_jmt(index3(1,1):index3(1,2),:)] * beta ...
                                                     - alpha .* p(index3(1,1):index3(1,2),:) + xi_jmt(index3(1,1):index3(1,2),:)...
                                                     + repmat(sigma(1),J-2,1) - p(index3(1,1):index3(1,2),:)*sigma(2); 
                                                                     % Since this
                                                                     % specific M*T
                                                                     % market lacks
                                                                     % one good
                                                                     % compared to
                                                                     % the other markets
                                                                      
for k = 2:M*T
                RC_exp(index2(k,1):index2(k,2),:)  = [ones(J-1,1) x_jmt(index3(k,1):index3(k,2),:)] * beta ...
                                                     - alpha .* p(index3(k,1):index3(k,2),:) + xi_jmt(index3(k,1):index3(k,2),:)...
                                                     + repmat(sigma(1),J-1,1) - p(index3(k,1):index3(k,2),:)*sigma(2);   
end

nom               = exp(RC_exp);
cumsum_delta_jmt  = cumsum(nom);
cumsum_1          = cumsum_delta_jmt(J-1,:);
totalsum          = diff(cumsum_delta_jmt(index,:)); 
totalsum          = [cumsum_1;totalsum];                                   % M*T total share of each j \in J_{mt}
Totalsum          = repelem(totalsum,J);                                   % Stretching it out to (J*M*T)*1 (since we have J-1 products in each market)
Totalsum          = Totalsum(2:end,:);

denom             = 1 + Totalsum;
s__jmt            = nom./denom;                                            % Getting ALL shares, i.e. including s0_mt 

inside_s_jmt      = zeros(Total-M*T-1,1);                                  % Now, getting shares of inside goods

for i = 1: M*T
inside_s_jmt(index3(i,1):index3(i,2))  = s__jmt(index2(i,1):index2(i,2),1);% Getting shares of inside goods
end


end
