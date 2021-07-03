function [index, indexx, index2 ,index3 ,m_index, m_indexx, m_index2, m_index3, m_drop_index] = indexes(J,M,T,Total,adj)

index   = [J:J:Total]';             % indexing each good in every market M*T combination
indexx  = [J-1:J-1:Total-M*T]';     % indexing that will be used in calculating outside good share in this setting

index2             = zeros(M*T,2);  % indexing the inside goods in each market in matrix form
index2(1,1)        = 1;
index2(1,2)        = J-1;

for i = 2 : M*T
    index2(i,1)    = index(i-1,1)+1;       
    index2(i,2)    = index(i,1)-1; 
end

index3             = zeros(M*T,2);  % indexing optimizer prices and shares etc.
index3(1,1)        = 1;
index3(1,2)        = J-1;

for i = 2 : M*T
   index3(i,1)     = index3(i-1,2)+1;
   index3(i,2)     = index3(i,1)+J-2;
end

m_index   = [J-adj:J-adj:Total-adj*M*T]';             % indexing each good in every market M*T combination for manipulated case

m_indexx  = [J-1-adj:J-1-adj:Total-M*T-adj*M*T]';     % indexing that will be used in calculating outside good share for manipulated case

m_index2             = zeros(M*T,2);                  % indexing the inside goods in each market in matrix form for manipulated case
m_index2(1,1)        = 1;
m_index2(1,2)        = J-1-adj;

for i = 2 : M*T
    m_index2(i,1)    = m_index(i-1,1)+1;       
    m_index2(i,2)    = m_index(i,1)-1; 
end

m_index3             = zeros(M*T,2);                  % indexing optimizer prices for manipulated case
m_index3(1,1)        = 1;
m_index3(1,2)        = J-adj-1;

for i = 2 : M*T
   m_index3(i,1)     = m_index3(i-1,2)+1;
   m_index3(i,2)     = m_index3(i,1) + J-adj-2;
end

mm_drop_index = zeros(adj,M);     % these elements will be dropped in manipulation: the dropped brands will be same for all time periods for each market

for i = 1 : M
mm_drop_index(:,i) = randperm(J-1,adj)';
end

mm_drop_index = [16,18,6,13,4,14,15,16,11,12,5,13,12,10,9,9,7,2,12,15;12,16,12,16,13,8,1,11,10,3,18,5,2,11,2,7,9,5,17,7;2,6,13,3,14,9,10,10,1,5,1,16,15,19,3,5,18,1,16,6;6,17,18,1,16,6,14,8,18,15,8,15,7,8,15,14,8,16,11,13;19,5,10,17,12,17,2,1,8,14,13,12,19,1,11,6,5,3,18,11];
mm_drop_index = sort(mm_drop_index);

m_drop_index  = zeros(adj,M*T);

t_index        = zeros(M,2);
t_index(1,1)   = 1;
t_index(1,2)   = T;

for  i = 2:M
   t_index(i,1) = t_index(i-1,2) + 1;
   t_index(i,2) = t_index(i,1) + T - 1;
end

for i = 1 : M
    m_drop_index(:,i:M:M*T-M+i) = repmat(mm_drop_index(:,i),1,T);
end





end
