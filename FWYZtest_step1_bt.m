% R is p by T return matrix, p is dimension and T is sample size
% K is the number of common factors
% Kernel_set is pre-calculating kernel weighting matrix
function [S_MT_set,S_MT,res_set] = FWYZtest_step1_bt(R,K,Kernel_set)
    RR = R'*R;
    [p,T] = size(R);
    [eigvector,eigvalue] = eig(RR);
    [~,dec_index] = sort(diag(eigvalue),'descend');
    eigK = eigvector(:,dec_index(1:K));
    Fhat = sqrt(T)*eigK';
    B = (1/T).*R*Fhat';
    res_set = R - B*Fhat;
    h = (2.35/sqrt(12))*T^(-0.2)*p^(-0.1);
    [S_MT] = FWYZtest_step1(R,K,Kernel_set);
    sigma_0  = res_set*res_set'./T;
    % conduct the bootstrap procedure for step 1 loading test
    for i = 1 : p
        for j = i+1 : p
           sigma_0(i,j) = sigma_0(i,j)*((1-0.01)^(j-i)); 
           sigma_0(j,i) = sigma_0(i,j);
        end
    end
    % the loop for bootstrap procedure is set to be 199
    S_MT_set = zeros(199,1);
    parfor Boot = 1 :199
        Y_star = B*Fhat + sqrtm(sigma_0)* mvnrnd(zeros(p,1),eye(p),T)';
        S_MT_set(Boot) = FWYZtest_step1(Y_star,K,Kernel_set);
    end
end