function [modularity,entropy] = evaluate(adj_mat,attr_mat,partition)
%
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
% 
% Description - evaluates the structure and attribute quality of the 
%               resulting clustering;
%               modularity is to measure structure quality and entropy is
%               to measure attribute quality for each attribute.
%
% Input  - adj_mat    : adjacency matrix
%        - attr_mat   : Tx1 cell array of attribute indicator matrices
%        - partition  : the final posterior distribution of the cluster labels
%
% Output - modularity : scalar
%        - entropy    : Tx1 vector of entropies
% -------------------------------------------------------------------------

    cltId = unique(partition);
    N = length(partition);
    K = length(cltId);   

	% -----------------------modularity------------------------
    M = sum(adj_mat(:));    
    
    rho = zeros(K);    
    clt = true(N,K);
    cltSize = zeros(K,1);
    for k = 1:K
        clt(:,k) = logical(partition==cltId(k));
        cltSize(k) = sum(clt(:,k));
        rho(k,k)   = sum(sum(adj_mat(clt(:,k),clt(:,k))))/M;
    end    
    for i = 1:(K-1)        
        for j = (i+1):K            
            rho(i,j) = sum(sum(adj_mat(clt(:,i),clt(:,j))))/M;
            rho(j,i) = rho(i,j);
        end
    end
    
    alpha = sum(rho,2);
    modularity = sum(diag(rho)-alpha.^2);
    
    % -------------------------entropy-------------------------
    N = length(partition);
    T = length(attr_mat);

    entropy = zeros(T,1);				
    coefficient = cltSize/N;			
    for t=1:T			
        for k=1:K
            distribution = sum(attr_mat{t}(clt(:,k),:))/cltSize(k);
            index = logical(distribution>0);
            ent   = -dot(distribution(index),log(distribution(index)))/log(2);
            entropy(t) = entropy(t) + ent * coefficient(k);
        end						
    end
end