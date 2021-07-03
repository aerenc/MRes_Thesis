function [m_x_jmt,m_w_jmt,m_xi_jmt,m_omega_jmt,m_mc_jmt]  = simulateDataManipulated(Total,J,M,T,adj,x_jmt,w_jmt,xi_jmt,omega_jmt,m_drop_index,gamma)

% Now, we drop 5 elements for each market for every time period randomly:

% First, we reshape the data into (J-1) * (M*T) matrices for each element of
% x and w and dropping the elements accordingly to m_drop_index that we
% constructed in INDEXES subsection above:

m_reshaped_x1_jmt    = reshape(x_jmt(:,1),J-1,M*T); 
m_reshaped_x2_jmt    = reshape(x_jmt(:,2),J-1,M*T); 
m_reshaped_x3_jmt    = reshape(x_jmt(:,3),J-1,M*T); 
m_reshaped_w1_jmt    = reshape(w_jmt(:,1),J-1,M*T); 
m_reshaped_w2_jmt    = reshape(w_jmt(:,2),J-1,M*T); 
m_reshaped_w3_jmt    = reshape(w_jmt(:,3),J-1,M*T); 
m_reshaped_xi_jmt    = reshape(xi_jmt,J-1,M*T); 
m_reshaped_omega_jmt = reshape(omega_jmt,J-1,M*T);

for  i = 1 : M*T
    m_reshaped_x1_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_x2_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_x3_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_w1_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_w2_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_w3_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_xi_jmt(m_drop_index(:,i),i) = 0;
end
for  i = 1 : M*T
    m_reshaped_omega_jmt(m_drop_index(:,i),i) = 0;
end

% m_ <=> manipulated:

m_x1_jmt = m_reshaped_x1_jmt(:);
m_x1_jmt = m_x1_jmt(m_x1_jmt~=0);
m_x2_jmt = m_reshaped_x2_jmt(:);
m_x2_jmt = m_x2_jmt(m_x2_jmt~=0);
m_x3_jmt = m_reshaped_x3_jmt(:);
m_x3_jmt = m_x3_jmt(m_x3_jmt~=0);

m_x_jmt = [m_x1_jmt m_x2_jmt m_x3_jmt] ;

m_w1_jmt = m_reshaped_w1_jmt(:);
m_w1_jmt = m_w1_jmt(m_w1_jmt~=0);
m_w2_jmt = m_reshaped_w2_jmt(:);
m_w2_jmt = m_w2_jmt(m_w2_jmt~=0);
m_w3_jmt = m_reshaped_w3_jmt(:);
m_w3_jmt = m_w3_jmt(m_w3_jmt~=0);

m_w_jmt = [m_w1_jmt m_w2_jmt m_w3_jmt];

m_xi_jmt = m_reshaped_xi_jmt(:);
m_xi_jmt = m_xi_jmt(m_xi_jmt~=0);

m_omega_jmt = m_reshaped_omega_jmt(:);
m_omega_jmt = m_omega_jmt(m_omega_jmt~=0);

m_mc_jmt = [ones(Total-M*T-adj*M*T,1) m_w_jmt] * gamma + m_omega_jmt;                    % getting marginal costs for manipulated data

end