function posterior=init(adj_mat,attr_mat,clt_num)
%
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
%
% Description - initializes the posterior distribution of the cluster
% labels Z by the same method as in Inc-Cluster.
%
% Input  - adj_mat   : 	adjacency matrix
%   	 - attr_mat  : 	Tx1 cell array of attribute indicator matrices												
%  		 - clt_num   : 	cluter number (K)
%
% Ouput  - posterior :  clustering of vertices by a simple method can be 
%                       taken as a rough approximation of posterior of Z;
%						the same initial clustering as in Inc-Cluster is 
%                       used for the fairness of comparsion
% -------------------------------------------------------------------------
   
    % --------------construct attribute augmented graph-------------
    N = size(adj_mat,1);
    T = length(attr_mat);
        
    vec = zeros(T+1,1);    
    vec(1) = N;    
    for t=1:T        
    	vec(t+1) = vec(t) + size(attr_mat{t},2);
    end
    
    n = nnz(adj_mat) + 2*T*N;
    aug_mat = spalloc(vec(end),vec(end),n);
    
    aug_mat(1:N,1:N) = adj_mat;
    for t=1:T
    	aug_mat(1:N,(vec(t)+1):vec(t+1)) = attr_mat{t};
    	aug_mat((vec(t)+1):vec(t+1),1:N) = attr_mat{t}';
    end
    
%     m = nnz(adj_mat);
%     n = m + 2*T*N;    
%     coo = zeros(n+1,3);
%     [coo(1:m,1),coo(1:m,2),coo(1:m,3)] = find(adj_mat);    
%     tuple = zeros(N,3);
%     s = m;
%     for t=1:T
%         [tuple(:,1),tuple(:,2),tuple(:,3)] = find(attr_mat{t});
%         coo((s+1):(s+N),1) = tuple(:,1);
%         coo((s+1):(s+N),2) = tuple(:,2)+vec(t);
%         coo((s+1):(s+N),3) = tuple(:,3);
%         s = s+N;        
%         coo((s+1):(s+N),1) = coo((s-N+1):s,2);
%         coo((s+1):(s+N),2) = coo((s-N+1):s,1);
%         coo((s+1):(s+N),3) = coo((s-N+1):s,3);
%         s = s+N;        
%     end
%     coo(end,:) = [vec(end),vec(end),0];
%     aug_mat = spconvert(coo);
%     clear coo tuple;
    
    % ------------------calculate transition matrix-----------------
    d = sum(aug_mat,2);
    d = spdiags(1./d,0,vec(end),vec(end));
    transition = d*aug_mat;
    clear aug_mat d;    

    % ------------------internal control parameters-----------------
    prune_size = 0.005;
    restart = 0.15;
    iter = 10;
    
    % ------------------neighborhood random walk distance-----------
    factor = (1-restart)*transition;
    clear transition;

    iterator = restart*factor;
    random_walk = iterator;    
    for i = 2:iter
        iterator = iterator*factor;
        value = logical(iterator>0.0 & iterator<=prune_size);
        iterator(value) = 0;
        clear value;
        iterator = sparse(iterator);
        random_walk = random_walk + iterator;
    end;
    clear iterator factor;
    
    % ---------------------centroids--------------------------------
    density = sum(adj_mat,2);										
    [~,index] = sortrows(density);
    center_label = index(N+1-(1:clt_num));
    clear density index;
    
    % ---------------------partition--------------------------------
    center_mat = random_walk(1:N,center_label);
    clear rand_walk;
    [~, partition] = max(center_mat,[],2);
    clear center_mat;
    partition(center_label) = 1:clt_num;
   
    % ---------------------posterior--------------------------------   
    posterior = sparse(1:N,partition,1,N,clt_num,N);
end