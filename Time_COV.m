% This code is used for estimating Time-varying covariance matrix proposed by us.
% R is p*T return matrix, p is dimension, T is sample size, 
% PCV is the sparse penalty coefficience for residual sparse matrix, 
% f_number is the number of factor
% tau is positive definite tunning parameter
% h is bandwidth parameter, default is rule of thumb

function [Sigma_r,Sigma_e,Residuals] = Time_COV(R,PCV,f_number,tau,h)
    p = size(R,1);
    T = size(R,2);
    if nargin == 4
        h = (2.35/sqrt(12))*(T^-0.2)*(p^-0.1);
    end
	Kernel_set = zeros(T,T);
	for  r = 1 :T
		for t = 1 :T
		   k = (h^-1) * Kernel2(T,t,r,h);
		   Kernel_set(t,r) = k^0.5; 
		end
    end
    % estimate time varying factor loading and the corresponding factors
	[Factor_loadings,Factor] = Timevarying_factor_model(R,f_number,Kernel_set,h);
 
    % obtain the residual
    Residuals = [];
    for t = 1 : T
        Factor_loading_t = Factor_loadings{t};
        Residual_t = R(:,t)- Factor_loading_t'*Factor(:,t);
        Residuals = [Residuals, Residual_t];
    end
    
    % save the residual, and perform the sparse estiamtor by using packages 'spcov' in R. 
    subdata = Residuals';
    P = PCV*(ones(p) - eye(p));
    tau = tau;
    save('Mat_R_temp.mat','subdata','P','tau');
    [s,e] = dos('Rspcov.bat'); 
    Sigma_e  = csvread('test.csv');
    Sigma_r = Factor_loading_t'* (cov(Factor'))*Factor_loading_t + Sigma_e; % our time-varying covariance matrix estimator
    mu = Factor_loading_t'* mean(Factor')';
end