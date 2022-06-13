% This program is used to test H0: NO time-varying factor loading H1 time-varying factor loading
% R is d*T return matrix, d is dimension, T is sample size, 
% p is the number of factor
% Kernel_set is pre-calculating kernel weighting matrix to improve the calculating speed
% For test procedure, we set h = 2.35/sqrt(12) T^{-0.2}N^{-0.1}
function [J_NT_hat,E_hat] = FWYZtest_step1(R,p,Kernel_set)
    %%%% calculate conventional PCA
    RR = R'*R;
    [N,T] = size(R);
    [eigvector,eigvalue] = eig(RR);
    [~,dec_index] = sort(diag(eigvalue),'descend');
    eigK = eigvector(:,dec_index(1:p));
    Ftilde = sqrt(T)*eigK';
    B = (1/T).*R*Ftilde';
    F = Ftilde;
    %%%
    h = (2.35/sqrt(12))*T^(-0.2)*N^(-0.1);
    [loading_hat,F_hat] = Timevarying_factor_model(R,p,Kernel_set,h);

    l_set = [];
    for j= 1:20
        l_temp = loading_hat{j};
        l_set(:,j) = l_temp(:,2);
    end
    
    M_hat = 0;
    for i = 1 :N
        for t =1 :T
            loading_hat_temp = loading_hat{t};
            M_hat = M_hat +  (loading_hat_temp(:,i)'*F_hat(:,t) - B(i,:)*F(:,t))^2;
        end
    end
    M_hat = M_hat/(N*T);
    
    % calulate residual
    E_hat = [];
    for t = 1: T
        Factor_loading_t =loading_hat{t};
        Residual_t = R(:,t)- Factor_loading_t'*F_hat(:,t);
        E_hat = [E_hat, Residual_t];
    end
    h1  = h;
    Kernel_set2 = Kernel_set.^2;
    
    Fs = F_hat';
    Fss = F';
    B_nt_hat_set = [];
    for i =1 :N
        Kcellall =  mat2cell(Kernel_set2,ones(T,1),ones(T,1));  
        Fscell =mat2cell(Fs,ones(T,1),p);
        Fscellall =   repmat(Fscell ,1,T);
        Ftcellall =  mat2cell(repmat(F_hat,T,1),p.*ones(1,T),ones(1,T));
        Fsscell =   mat2cell(Fss,ones(T,1),p);
        Fsscellall =   repmat(Fsscell ,1,T);
        Fttcellall =  mat2cell(repmat(F,T,1),p.*ones(1,T),ones(1,T));
        ei2 = num2cell((E_hat(i,:).^2)');
        ei2all = repmat(ei2,1,T);
        KxFsall = cellfun(@mtimes,Kcellall,Fscellall,'UniformOutput',false);
        KxFsxFtall = cellfun(@mtimes,KxFsall,Ftcellall,'UniformOutput',false);             
        FsFtall = cellfun(@mtimes,Fsscellall,Fttcellall,'UniformOutput',false);
        brakall =  cellfun(@minus,KxFsxFtall,FsFtall,'UniformOutput',false);
        brak2all =  cellfun(@times,brakall,brakall,'UniformOutput',false);
        partall = cellfun(@mtimes,brak2all,ei2all,'UniformOutput',false);
        B_nt_hat_set(i) =sum(sum(cell2mat(partall)));
    end
    B_nt_hat = sum(B_nt_hat_set);
    B_nt_hat = B_nt_hat*(h1^0.5)/((T^2) * (N^0.5));

   
    Sigma_F_hat = 1/T *F_hat *F_hat';
    V_nt_hat = 0;
    for s = 1 :T
        for r = 1:T
            if s == r
                continue;
            end
            V_nt_hat = V_nt_hat + kernel_bar(T,s,r,h1)^2 .* ((F_hat(:,s)'*Sigma_F_hat*F_hat(:,r))^2)*(E_hat(:,r)'*E_hat(:,s))^2;
        end
    end
    V_nt_hat = 2* T^(-2)*N^(-1)*h1^(-1)*V_nt_hat;

    J_NT_hat = V_nt_hat^(-0.5)*(T*(N^(0.5)) *(h1^0.5)*M_hat -B_nt_hat);

end