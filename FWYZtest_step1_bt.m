function [S_MT_set,S_MT,res_set] = FWYZtest_step1_bt(Y,K,Kernel_set)
    RR = Y'*Y;
    [N,T] = size(Y);
    [eigvector,eigvalue] = eig(RR);
    [~,dec_index] = sort(diag(eigvalue),'descend');
    eigK = eigvector(:,dec_index(1:K));
    Fhat = sqrt(T)*eigK';
    B = (1/T).*Y*Fhat';
    res_set = Y - B*Fhat;
    h = (2.35/sqrt(12))*T^(-0.2)*N^(-0.1);
    [S_MT] = FWYZtest_step1(Y,K,Kernel_set);
    sigma_0  = res_set*res_set'./T;
    % conduct the bootstrap procedure for step 1 loading test
    for i = 1 : N
        for j = i+1 : N
           sigma_0(i,j) = sigma_0(i,j)*((1-0.01)^(j-i)); 
           sigma_0(j,i) = sigma_0(i,j);
        end
    end
    % the loop for bootstrap procedure is set to be 199
    S_MT_set = zeros(199,1);
    parfor Boot = 1 :199
        Y_star = B*Fhat + sqrtm(sigma_0)* mvnrnd(zeros(N,1),eye(N),T)';
        S_MT_set(Boot) = FWYZtest_step1(Y_star,K,Kernel_set);
    end
end