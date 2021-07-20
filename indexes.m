function [base_index, base_indexx, base_index2 ,base_index3, index, indexx, index2, index3, indexxx, IDmkt] = indexes(J,M,T,Total)

base_index   = [J:J:Total]';                                               % Indexing each good in every market M*T combination
base_indexx  = [J-1:J-1:Total-M*T]';                                       % Indexing that will be used in calculating outside good share in this setting

base_index2             = zeros(M*T,2);                                    % Indexing the inside goods in each market in matrix form
base_index2(1,1)        = 1;
base_index2(1,2)        = J-1;

for i = 2 : M*T
    base_index2(i,1)    = base_index(i-1,1)+1;       
    base_index2(i,2)    = base_index(i,1)-1; 
end

base_index3             = zeros(M*T,2);                                    % Indexing optimizer prices and shares etc.
base_index3(1,1)        = 1;
base_index3(1,2)        = J-1;

for i = 2 : M*T
   base_index3(i,1)     = base_index3(i-1,2)+1;
   base_index3(i,2)     = base_index3(i,1)+J-2;
end


%% "REAL" DATA INDEXES:

index  = [J-1:J:J*M*T-1]';                                                 % Indexing each good in every market M*T combination for "real" data
indexx = [J-2:J-1:(J-1)*M*T-1]';                                           % Indexing that will be used in calculating outside good share in this setting for "real" data


index2 = zeros(M*T,2);                                                     % Indexing the inside goods in each market in matrix form for "real" data
index2(1,1) = 1;
index2(1,2) = J-2;
for k = 2 : M*T
    index2(k,1) =  (k-1)*J;
    index2(k,2) =  (k-1)*J+J-2;
end

index3 = zeros(M*T,2);                                                     % Indexing optimizer prices and shares etc. for "real" data
index3(1,1) = 1;
index3(1,2) = J-2;
for k = 2 : M*T
             index3(k,1) = (J-1)*(k-1);
             index3(k,2) = (J-1)*(k-1)+J-2; 
end


% Now produce a matrix which indicates the number of products in each
% market (<=> "prods" in EIO1 assignment):

indexxx = zeros(M*T,1);
indexxx(1,:) = indexx(1,:) ;

for   k = 2 : M*T
   indexxx(k,:) = indexx(k,:) - indexx(k-1,:);
end


IDmkt = [1:1:M*T]';

IDmkt = repelem(IDmkt,J-1);

IDmkt = IDmkt(2:end,:);






end
