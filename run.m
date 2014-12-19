function run(adj_coo_file, attr_tab_file, clt_num, iter_num, init_clt_file)
%
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
%
% Description - provide a single command for the use of BAGC
%
% Input  - adj_coo_file  : file name for adjacency matrix in format COO 
%                          (coordinate list)
%        - attr_tab_file : file name for attribute table with a row for 
%                          each vertex and a column for each attribute
%        - clt_num       : cluster number
%        - iter_num      : iteration number, optional
%        - init_clt_file : file name for initial clustering (partition, or 
%                          assignment) of vertices;
%                          vector of size N after loaded in which the number
%                          of distinct elements should be less than or equal
%                          to clt_num;
%                          optional argument.
% ------------------------------------------------------------------------

    % load data
    adj_coo  = load(adj_coo_file);
    attr_tab = load(attr_tab_file);
    
	% run bagc
    if nargin == 5
        [modularity,entropy,time] = bagc(adj_coo,attr_tab,clt_num,iter_num,load(init_clt_file));
    elseif nargin == 4
        [modularity,entropy,time] = bagc(adj_coo,attr_tab,clt_num,iter_num);
    elseif nargin == 3
        [modularity,entropy,time] = bagc(adj_coo,attr_tab,clt_num);
    else
        error('incorrect #input arguments!');
    end
    
    % display results
    disp('modularity:');            disp(modularity);
    disp('entropy:');               disp(entropy');
    disp('optimizing time(s):');    disp(time);
end