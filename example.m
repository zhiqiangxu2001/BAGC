function example()
%
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
% 
% Description - provide an example to illustrate the use of the code on the
%               given attributed graph. The output clustering has a 
%               modularity value 0.5513, and entropies 0.2318 for attribute
%               1 (column 1 of attr_tab.txt) and 0.3849 for attribute 2 
%               (column 2 of attr_tab.txt). The running time is about 
%               0.005s on a machine with 32-bit windows OS, 2.67GHz 4-core 
%               CPU and 4G memory.
%
%
% Output - modularity : scalar, structure quality measure
%        - entropy    : Tx1 entropy vector of entropies with an element 
%                       for each attribute; attribute quality measure.
%        - time       : optimizing time
% ------------------------------------------------------------------------

    % load data
    adj_coo  = load('adj_coo.txt');
    attr_tab = load('attr_tab.txt');

    % set cluster number    
    clt_num  = 3; 
    
    % set iteration number for optimization
    iter_num = 5;
    
	% run bagc
    [modularity,entropy,time] = bagc(adj_coo,attr_tab,clt_num,iter_num);
    
    disp('modularity:');            disp(modularity);
    disp('entropy:');               disp(entropy');
    disp('optimizing time(s):');    disp(time);

end