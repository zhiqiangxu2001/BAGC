function [adj_mat,attr_mat] = formatted(adj_coo,attr_tab)
%
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
%
% Description - convert raw data, i.e., coordinate list of adjacency matrix
% and attribute table, into internal representations.
%
% Input  - adj_coo  :  adjacency matrix in format COO (coordinate list);
%				 	   note that Matlab array indices are "1-based";
%                      add an additional row [N,N,0] to the list to ensure 
%                      a square matrix;
%                      other formats are also feasible only if they can be 
%                      easily converted into the final format, i.e., 
%                      adjacency matrix, after a bit of modification to the
%                      current code.
%
%		- attr_tab	:  NxT attribute table with a row for each vertex and a
%                      column for each attribute;
%				 	   other formats are also feasible only if they can be 
%                      easily converted into the final format shown below 
%                      after a bit of modification to the current code.
%
% Output - adj_mat  :  adjacency matrix in sparse representation
%
%		 - attr_mat :  Tx1 cell array of attribute indicator matrices;
%					   each attribute a^t corresponds to an NxM^t indicator
%                      matrix Y^t in which each vertex corresponds to a row,
%                      i.e., a one-of-M^t indicator vector.
% -------------------------------------------------------------------------

	% convert coordinate list to sparse matrix for adjacency matrix
	adj_mat  = spconvert(adj_coo);
		
    % convert attribute table to Tx1 cell array of attribute indicator
    % matrices
	N = size(adj_mat,1);
	T = size(attr_tab,2);
	
	attr_mat = cell(T,1);
    for t=1:T
        [V,~,J] = unique(attr_tab(:,t));
        attr_mat{t} = sparse(1:N,J,1,N,length(V),N);
    end
end