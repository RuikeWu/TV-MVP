% Solving mininmum variance problem
function [w] = Markowitz_MVP(R_cov)
    p = size(R_cov,1);
    N1 = ones(p,1);
    cvx_begin quiet
         variable w(p)
         minimize( w'*R_cov*w)
         subject to
             N1' *w == 1 
     cvx_end
end