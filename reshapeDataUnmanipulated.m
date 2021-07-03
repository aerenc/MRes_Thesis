function [array_s_jmt,array_p_jmt,array_xi_jmt,array_omega_jmt,array_mc_jmt,array_x_jmt,array_w_jmt,reshaped_x1_jmt,reshaped_x2_jmt,reshaped_x3_jmt,...
    reshaped_w1_jmt,reshaped_w2_jmt,reshaped_w3_jmt] = reshapeDataUnmanipulated(J,M,T,s_jmt,p_jmt,xi_jmt,omega_jmt,mc_jmt,x_jmt,w_jmt)


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


reshaped_x1_jmt    = reshape(x_jmt(:,1),J-1,M*T);
reshaped_x2_jmt    = reshape(x_jmt(:,2),J-1,M*T);
reshaped_x3_jmt    = reshape(x_jmt(:,3),J-1,M*T); 
reshaped_w1_jmt    = reshape(w_jmt(:,1),J-1,M*T); 
reshaped_w2_jmt    = reshape(w_jmt(:,2),J-1,M*T); 
reshaped_w3_jmt    = reshape(w_jmt(:,3),J-1,M*T); 
reshaped_xi_jmt    = reshape(xi_jmt,J-1,M*T); 
reshaped_omega_jmt = reshape(omega_jmt,J-1,M*T);

end