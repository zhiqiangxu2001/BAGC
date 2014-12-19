function [modularity,entropy,time] = bagc(adj_coo,attr_tab,clt_num,iter_num,init_clt)
% 
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
%
% Description - main function to implement the proposed algorithm BAGC.
%
% Input  - adj_coo    : adjacency matrix in format COO (coordinate list)
%        - attr_tab   : attribute table with a row for each vertex and a
%                       column for each attribute
%        - clt_num    : cluster number
%        - iter_num   : iteration number, default 5
%        - init_clt   : vector of size N in which the number of distinct 
%                       elements is less than or equal to clt_num;
%                       an optional argument as initial clustering (partition,
%                       or assignment) of vertices privided by users.
% 
% Output - modularity : scalar, structure quality measure
%        - entropy    : Tx1 entropy vector with an element for each attribute;
%                       attribute quality measure.
%        - time       : optimizing time
% -------------------------------------------------------------------------

	% format data --------------------------------
	[adj_mat,attr_mat] = formatted(adj_coo,attr_tab);
	clear adj_coo attr_tab;
	
	% init ---------------------------------------
    if nargin == 5
        N = length(init_clt);
        [V,~,J] = unique(init_clt);
        assert(length(V)<=clt_num,'#labels of clusters exceeds specified #clusters');
        posterior = sparse(1:N,J,1,N,clt_num,N);
    elseif nargin == 4
        posterior = init(adj_mat,attr_mat,clt_num);
    elseif nargin == 3
        posterior = init(adj_mat,attr_mat,clt_num);
        iter_num  = 5;
    else
        error('incorrect #input arguments!');
    end
	
	% optimize -----------------------------------
	[posterior,time] = optimize(adj_mat,attr_mat,clt_num,iter_num,posterior);
    
    % partition - clustering of vertices ---------
	[~,partition] = max(posterior,[],2);   
    
	% evaluate -----------------------------------
	[modularity,entropy] = evaluate(adj_mat,attr_mat,partition);
end