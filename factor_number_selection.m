% R is the input return matrix , N X T matrix, N is dimension, T is sample size
% max_number is the maximum number of factor
% This code use the BIC-type information criterion to select the number of factors in appendix B.3
function [f_number,IC_R_set] = factor_number_selection(R,max_number)
    IC_R_set = [];
    [N,T] = size(R);
    h =(2.35/sqrt(12))*(T^-0.2)*(N^-0.1);
    Kernel_set = zeros(T,T);
	for  r = 1 :T
		for t = 1 :T
		   k = (h^-1) * Kernel2(T,t,r,h);
		   Kernel_set(t,r) = k^0.5; 
		end
	end
    
    for f_count = 1 :max_number
        [Factor_loadings,Factor,X_r_cell,~,h] = Timevarying_factor_model(R,f_count,Kernel_set,h);
        Lambda_wan_set = cell(T,1);
        for r = 1 :T
            X_r_temp = X_r_cell{r};
            Lambda_r_temp = Factor_loadings{r};
            Lambda_wan_set{r} = 1/(N*T) * X_r_temp'*X_r_temp*Lambda_r_temp';
        end
        m = size(Factor,1);
        cvx_begin quiet
            variable F(m,T)
            expression FF((m*T),T)
            minimize V(R,F,Lambda_wan_set,FF)
        cvx_end            
        V_R = (V(R,F,Lambda_wan_set,FF)^2)/(N*T);
        rou = (N+ T*h)/(N*T*h)*log((N*T*h)/(N+ T*h)); 
        IC_R = log(V_R) + rou* m; 
        IC_R_set = [IC_R_set;IC_R];
    end
    [~,f_number] = min(IC_R_set);
end