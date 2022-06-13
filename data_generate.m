function [Y,r_cov,s_cov]=data_generate(d,T,rho)
    rng('default');rng(2)
    subt = [1:T]/T; % generate scaled time

    %generate the smooth time varying factor loading defined in appendix C.2.1
    loading1=[];
    loading2 =[];
    for subi = 1:d
       z =  (rho).*((3+randn(1))*subt + randn(1).*sin(4*pi*subt))+randn(1);
       z2 =  (rho).*((3+randn(1))*subt + randn(1).*sin(4*pi*subt)+(randn(1)*subt).^2)+randn(1);
       loading1(subi,:) = z;
       loading2(subi,:) = z2;
    end

    %cross-sectional dependence in residuals
    s_cov = [];
    for i = 1 : d
    for j = i+1 : d
       s_cov(i,j) = ((0.5)^abs(j-i)); 
       s_cov(j,i) = s_cov(i,j);
    end
    end
    s_cov = s_cov +eye(d);                    
    r_cov = [loading1(:,T),loading2(:,T)]*[loading1(:,T),loading2(:,T)]'+s_cov;
    
    F1_set = [];
    F2_set = [];
    F1_t0 = randn(1);
    F2_t0 = randn(1);
    for j = 1:T+200
          F1_t = 0.6*F1_t0 + mvnrnd(0,0.64,1);
          F1_set(j) = F1_t;
          F1_t0 = F1_t;
          F2_t = 0.3*F2_t0 + mvnrnd(0,0.91,1);
          F2_set(j) = F2_t;
          F2_t0 = F2_t;
    end
    F =  [F1_set(201:200+T);F2_set(201:200+T)];

    x_k = mvnrnd(zeros(d,1),s_cov,T)' ; 

    fa_part = [];
    for qq = 1 : T
        B1_temp = [loading1(:,qq),loading2(:,qq)];
        fa_part(:,qq) = B1_temp*F(:,qq);
    end
    Y = fa_part+ x_k;
end