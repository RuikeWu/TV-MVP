% R is the input return matrix, dimension times sample size
% lam is sparse tuning parameter
% K is the factor number
% tau is positive definite tunning parameter
% This code selects the value for bandwidth by using cross-validation method in appendix B.2

function [h_cv,All_d] = CV_for_h(R,lam,K,tau)

    [P,T] = size(R);
    N1 = ones(P,1);
    All_d = []; 
    h = [0.05:0.025:0.95];
    if P/T + 0.1 > 0.7
        p_T  = P/T;
        test_loop = round((1 - P/T-0.1)/0.1);
    else
        p_T = 0.6;
        test_loop = 3;
    end
    testing_set_cell = cell(test_loop,1);
    for kk = 1 : test_loop
        if p_T+0.1*(kk+1)> 1
            testing_set_cell{kk} = R(:,floor(T*(p_T+0.1*(kk)))+1:floor(T));
        else
            testing_set_cell{kk} = R(:,floor(T*(p_T+0.1*(kk)))+1:floor(T*(p_T+0.1*(kk+1))));
        end
    end
    training_set_cell = cell(test_loop,1);
    for kk = 1 : test_loop
        training_set_cell{kk} = R(:,1:floor(T*(p_T+0.1*(kk))));
    end
    SR_temp = [];
    for i =1:length(h)
            for j = 1 : test_loop
                fprintf('CV for h Large loop %d Num. %d\n',i,j);
                training_set = training_set_cell{j};
                testing_set = testing_set_cell{j};
                [Sigma] = Time_COV(training_set,lam,K,tau,h(i));
                C = N1' * Sigma * N1;
                wopt = (1/C) .* (Sigma * N1);
                real_return = [];
                for subi = 1 :size(testing_set,2)
                   real_return = [real_return;wopt'* testing_set(:,subi)];
                end
                SR1 = mean(real_return)/std(real_return);
                SR_temp(i,j)= SR1;
            end
    end
     All_d = sum(SR_temp,2);
    [~,loc] = max(All_d);
    h_cv = h(loc);
end


