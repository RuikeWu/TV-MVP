% R_history is the input data, dimension times sample size
% k is the the folds
% tau is positive definite tunning parameter
% This code selects the value for tuning parameter lambda by using cross-validation method in appendix B.1
function [lambda,All_d] = CV_for_spcov(R_history,k,tau)

    randn('state',12345);rand('state',1234567); 
    P0 = [0.05:0.05:1.5,10:10:100,150:50:1000,1250:250:6000];
    All_d = []; 
    rand_rank = randperm(size(R_history',1)); %generate random number without replications
    subnum = floor(size(R_history,2)/k);
    for i =1:length(P0)
            P = P0(i);
            Diff = 0;
            P = P*(ones(size(R_history,1)) - eye(size(R_history,1)));
            tau = tau;
            for  j = 1 : k
                testing_data = R_history(:,rand_rank((j-1)*subnum+1:j*subnum));
                training_data = R_history;
                training_data(:,rand_rank((j-1)*subnum+1:j*subnum)) = [];
                subdata = training_data';
                save('Mat_R_temp.mat','subdata','P','tau');
                fprintf('Getting SPCOV under %.3f cross Num. %d\n',P0(i),j);
                [mm,mmm]=dos('Rspcov.bat');
                sp_cov = csvread('test.csv');
                S_v = cov(testing_data');
                if  ~isreal(log_det(sp_cov))   %log_det is a function in CVX toolbox
                    Dt = - log(det(sp_cov)) - trace(S_v/sp_cov);    
                else
                    Dt = - log_det(sp_cov) - trace(S_v/sp_cov);
                end
                Diff = Diff + Dt;
            end
            All_d(i) = Diff/ k;
            if length(nonzeros(sp_cov)) == size(R_history,1)
                break;
            end
    end
    [~,loc] = max(All_d);
    lambda = P0(loc);
end