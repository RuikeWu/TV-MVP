% p is dimension
% T is sample size
% loop is replications for test
% alpha is the significance level
% This code runs for the size of DGP2(heteroskedasticity in residual) in
% our setting, for other types of DGP, one only need to change the way of
% generating virtual data.
function  [level]= size_test(p,T,loop,alpha)

%%%%%%%%%  Preliminary setting
    fprintf('Preliminary setting\n');
    K =2 ;
    B1_t = [];
    B2_t = [];
    % constant loading
    for j = 1:T+200
      B1_t(j) = randn(1);
      B2_t(j) = randn(1);
    end
    B1 = [B1_t(201+T - p:200+T);B2_t(201+T - p:200+T)]';
    % heteroskedasticity
    sigma_dset= [];
    for mm = 1 :p
              sigma_d = rand(1) + 0.5;
              sigma_dset(mm) = sigma_d;
    end
%%%%%%%%%% Generate virtual data
    fprintf('Generate virtual data\n');
    Y_t_set = cell(loop,1); 
    for i = 1:loop   
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
          for mm = 1 :p
              for kk = 1 :T
                e_it(mm,kk)= randn(1) *sigma_dset(mm);
              end
          end
         xk_set{i} = e_it;
         Y_t_temp = (B1*F)+ e_it;
         Y_t_set{i} = Y_t_temp;
    end
%%%%%% setting for  test

    p_val_set1= [];
    L_n_set1 =[];
    Stat = [];
    S_MT_boot_cell = cell(500,1);
    E_t_cell = cell(500,1);
    h = (2.35/sqrt(12))*T^(-0.2)*p^(-0.1);
    Kernel_set = zeros(T,T);
    for  r = 1 :T
        for t = 1 :T
           k = (h^-1) * Kernel2(T,t,r,h);
           Kernel_set(t,r) = k^0.5;
        end
    end

    

    parfor i = 1 :loop
        i
        Y = Y_t_set{i};
        tic;
        [S_MT_set,S_MT,E_t_temp] = FWYZtest_step1_bt(Y,K,Kernel_set);
        toc;
        E_t_cell{i} = E_t_temp;
        Stat(i) = sum(S_MT < S_MT_set)/199;
        S_MT_boot_cell{i} = S_MT_set;
    end

    parfor i = 1 :loop
        if Stat(i) < alpha(1)/2
            p_val_set1(i) = 0;
            L_n_set1(i) = -99;
            continue;
        end
        xk = E_t_cell{i};
        [Ln,ifa]= FWYZtest_step2(xk,alpha(1));
        p_val_set1(i) = ifa;
        L_n_set1(i) = Ln;
    end

    sz1 = sum(Stat<alpha(1)/2)/loop;
    sz2 = sum(p_val_set1)/(loop-sum(Stat<alpha(1)/2));
    sum_all = 0;
    for i = 1 :length(Stat)
        if Stat(i) < alpha(1)/2
            sum_all = sum_all +1;
            continue;
        else
            sum_all = sum_all + p_val_set1(i);
        end  
    end
    sz3 = sum_all/loop;
    level = [sz1,sz2,sz3];
end