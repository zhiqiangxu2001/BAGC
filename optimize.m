function [posterior,time]=optimize(adj_mat,attr_mat,clt_num,iter_num,posterior)
%
% Author - Zhiqiang Xu, 05/2012
%
% Email  - zxu1@e.ntu.edu.sg
%
% Description - iteratively optimizes the posterior of the cluster labels Z.
%
% Input  - adj_mat   : 	adjacency matrix
%        - attr_mat  : 	Tx1 cell array of attribute indicator matrices												
%  		 - clt_num   : 	cluter number (K)
% 		 - iter_num  : 	maximal number of external loops
% 		 - posterior : 	initial posterior of cluster labels Z;
%						NxK nonnegative matrix with row sum 1 assigned to 
%                       tilde_beta
% Output - posterior :  final posterior of cluster labels Z
%		 - time      : 	optimizing time
%--------------------------------------------------------------------------

 	% ----------internal control parameters----------
    log_min  = log(realmin);	
    iter_ext = iter_num;
	eps_ext  = 1e-10;       																						% stopping precision for external loop
    iter_int = 5;																									% maximal number of internal loop iterations
    eps_int  = 1e-4;                                                                                                % stopping precision for internal loop               
    prune    = 1e-10;           					                        										% prune size for tilde_beta used to speedup                   

    % ----------common constants---------------------
	N = size(adj_mat,1);
	T = length(attr_mat);
	K = clt_num;

	% ----------hyperparameters----------------------
    xi    = ones(1,K);               																				% hyperparameter of p(\alpha|\xi)
    mu	  = ones(K);               																					% hyperparameter of p(\phi|\mu,\nu)
    nu	  = ones(K);																								% hyperparameter of p(\phi|\mu,\nu)
    gamma = cell(T,1);										
    for t=1:T
        gamma{t} = ones(K,size(attr_mat{t},2));                                                                     % hyperparameter of p(\theta|\gamma)
    end        
                                
    % ----------optimize-----------------------------
    tic;  
              
    c = gammaln(sum(xi))-sum(gammaln(xi));                                                                          % constant in lower bound
    b = betaln(mu,nu);
    c = c-(sum(diag(b))+sum(b(:)))/2;                
    for t=1:T
        c = c+sum(gammaln(sum(gamma{t},2)))-sum(sum(gammaln(gamma{t})));
    end
    
    tilde_gamma = cell(T,1);
    tilde_beta  = posterior;																						% parameters of approximate distributions q(Z_i|\tilde{\beta}), i=1,...,N     
    
    i = 1; 
    while i < iter_ext                                                                                              % external loop       
        % update tilde_xi, tilde_mu, tilde_nu,tilde_gamma
        col_sum  = sum(tilde_beta);
        tilde_xi = xi + col_sum;																				    % hyperparameter of q(\alpha|\tilde{\xi})
        tilde_mu = tilde_beta'*adj_mat*tilde_beta;                                                                  % hyperparameter of q(\phi|\tilde{\mu},\tilde{\nu})
        tilde_nu = col_sum'*col_sum-tilde_beta'*tilde_beta-tilde_mu;                                                % hyperparameter of q(\phi|\tilde{\mu},\tilde{\nu})
        tilde_mu((K+1)*(1:K)-K) = tilde_mu((K+1)*(1:K)-K)/2;
        tilde_nu((K+1)*(1:K)-K) = tilde_nu((K+1)*(1:K)-K)/2;                        
        tilde_mu = mu + tilde_mu;
        tilde_nu = nu + tilde_nu;            
        for t=1:T
            tilde_gamma{t} = gamma{t}+tilde_beta'*attr_mat{t};                                                      % hyperparameter of p(\theta|\gamma)
        end
        % compute the lower bound
        d = -gammaln(sum(tilde_xi))+sum(gammaln(tilde_xi));
        b = betaln(tilde_mu,tilde_nu);
        d = d+(sum(diag(b))+sum(b(:)))/2;            
        for t=1:T
            d = d+sum(sum(gammaln(tilde_gamma{t})))-sum(gammaln(sum(tilde_gamma{t},2)));
        end
        d = d-sum(tilde_beta(logical(tilde_beta>0)).*log(tilde_beta(logical(tilde_beta>0)))); 
        if i > 1       
            old_lower_bound = lower_bound; 
        end
        lower_bound = c+d;																				        	% lower bound of the mariginal likelihood, i.e., log(P(X,Y))
        if i > 1
            if old_lower_bound==0
                error('divisor is 0');
            end
            err = abs(1-lower_bound/old_lower_bound);                                                               % relative error of lower bound
            if err < eps_ext 
                break;
            end
        end        
        % update tilde_beta
        b  = psi(tilde_nu);
        B1 = psi(tilde_mu)-b;
        B2 = b-psi(tilde_mu+tilde_nu);                        
        A  = repmat(psi(tilde_xi)-psi(sum(tilde_xi)),N,1);
        for t=1:T
            a = psi(tilde_gamma{t})-repmat(psi(sum(tilde_gamma{t},2)),1,size(attr_mat{t},2));
            A = A + attr_mat{t}*a';
        end
        diff = eps_int+1;   																						% absolute error of tilde_beta
        keep_on  = true;          
        j = 1;                       
        while (diff > eps_int) && (j < iter_int) && keep_on                                                         % internal loop
            old_tilde_beta = tilde_beta;                                   
            tmp = tilde_beta*B2';     
            tmp = repmat(sum(tmp),N,1)-tmp;
            tilde_beta = adj_mat*tilde_beta*B1'+tmp+A;                               
            tilde_beta(tilde_beta<log_min) = log_min;                
            tilde_beta = tilde_beta-repmat(max(tilde_beta,[],2),1,K);                                               % subtract the row maximum from each row
            tilde_beta = exp(tilde_beta);                     
            tilde_beta = tilde_beta./repmat(sum(tilde_beta,2),1,K);                                                 % normalize each row                
            tilde_beta(tilde_beta<prune) = 0;                
            if j > 1
                old_diff = diff;
                diff = max(max(abs(tilde_beta-old_tilde_beta)));                
                if diff >= old_diff
                    keep_on = false;                                                                                % diff increase, stop iterating
                end
            elseif j == 1                
                diff = max(max(abs(tilde_beta-old_tilde_beta)));                        
            end
            j = j+1;
        end
        i = i+1;
    end    
    time = toc;

    % ----------result-------------------------------
    posterior = tilde_beta;       
end