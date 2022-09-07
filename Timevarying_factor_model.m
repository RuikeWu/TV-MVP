% This function excute the local principal component method to simultaneously estimate 
% time varying factor loading and corresponding factors. 
% R is p X T matrix of excess return, p is dimension, T is sample size 
% f-number is the number of factor
% Kernel_set is pre-calculating kernel weighting matrix to improve the calculating speed
% h is the bandwidth used for Kernel_set calculating
function [Factor_loadings,Factor_set,X_r_cell,Fr_cell,h] = Timevarying_factor_model(R,f_number,Kernel_set,h)

    
    [p,T] = size(R);
    Factor_loadings = cell(T,1);
    Fr_cell = cell(T,1);
    X_r_cell = cell(T,1);
    
    for  r = 1:T
        Kernel_r = Kernel_set(:,r);
        X_r = [];
        for i = 1 :p
            X_i_r = Kernel_r .* (R(i,:)');
            X_r = [X_r,X_i_r];
        end
        X_r_cell{r} = X_r;
        XX_r = X_r * X_r';
        [XX_evec,XX_evalue] = eig(XX_r);
        [XX_evalue_dec,dec_index] = sort(diag(XX_evalue),'descend');
        %
        F_r = T^0.5 .*XX_evec(:,dec_index(1:f_number));
        
		% This long ¡®if¡¯ code regulizes the direction of eigenvector
         if r == 1
            rank1 =(abs(F_r(2,:))>abs(F_r(1,:)));
            two_rows = [F_r(1,:);F_r(2,:)];
            direction = [];
            for d_count = 1:length(rank1)
                direction = [direction,(two_rows(1,d_count)^(1-rank1(d_count)))*(two_rows(2,d_count)^rank1(d_count))];
            end
            direction = sign(direction);
        else
            if r <T
                two_rows = [F_r(r-1,:);F_r(r,:)];
                direction_else = [];
                for d_count = 1:length(rank1)
                    direction_else = [direction_else,(two_rows(1,d_count)^(1-rank1(d_count)))*(two_rows(2,d_count)^rank1(d_count))];
                end
                direction_else = sign(direction_else);


                for dir_c = 1:length(direction)
                     if direction_else(dir_c) ~= direction(dir_c)
                         F_r(:,dir_c) = -F_r(:,dir_c);
                     end
                end

                 rank1 =(abs(F_r(r+1,:))>abs(F_r(r,:)));
                two_rows = [F_r(r,:);F_r(r+1,:)];
                direction = [];
                for d_count = 1:length(rank1)
                    direction = [direction,(two_rows(1,d_count)^(1-rank1(d_count)))*(two_rows(2,d_count)^rank1(d_count))];
                end
                direction = sign(direction);
            else
                two_rows = [F_r(r-1,:);F_r(r,:)];
                direction_else = [];
                for d_count = 1:length(rank1)
                    direction_else = [direction_else,(two_rows(1,d_count)^(1-rank1(d_count)))*(two_rows(2,d_count)^rank1(d_count))];
                end
                direction_else = sign(direction_else);
                for dir_c = 1:length(direction)
                     if direction_else(dir_c) ~= direction(dir_c)
                         F_r(:,dir_c) = -F_r(:,dir_c);
                     end
                end
            end
        end
        %
        
        Fr_cell{r} = F_r;
        Lambda_r = 1/T .* (F_r'*X_r); 
        Factor_loadings{r} = Lambda_r;
    end

    
    Factor_set = [];
    for t= 1 :T
        F_t_hat_part1 = zeros(f_number,f_number);
        F_t_hat_part2 = zeros(f_number,1);
        factor_loading  = Factor_loadings{t};
        for i = 1 :p
            F_t_hat_part1 = F_t_hat_part1+ factor_loading(:,i)*factor_loading(:,i)';
            F_t_hat_part2 = F_t_hat_part2+ factor_loading(:,i)*R(i,t);
        end
        F_t_hat =  inv(F_t_hat_part1)*F_t_hat_part2;
        Factor_set =[Factor_set,F_t_hat];
    end
    
end