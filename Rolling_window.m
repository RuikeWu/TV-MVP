% This code achieve rolling window procedure in empirical analysis
% Data_matrix is the input matrix, the number of stocks by sample size
% history is the length of data you can get for history 
% persist is the time you hold on, e.g. weekly 5,  monthly 21
% P1 is tunning parameter for sparse residual covariance estimation
% cur_pos: star time point
% end_pos: end time point

% Output:  
% info: information set with Sharpe ratio, mean return and  standard deviation of propsoed strategy(loc 1,3,5) and '1/N'(loc 2,4,6)
% PD_return: realized portfolio return of TV-MVP
% eq_return: realized portfolio return of '1/N'
function [info,PD_return,eq_return] = Rolling_window(Data_matrix,history,persist,cur_pos,end_pos)
    
    tic;
    tau = 0.001;
    PD_return = [];
    eq_return = [];
    h_set = [];
    PCV_set = [];
    f_number_set =[];
    counth= 0;
    h = 0.5; 
    W_opt_set = cell(1);
    W_count = 1;
    try_time = 0;
    fprintf('Loop: ');
    for i = cur_pos:persist:end_pos-persist
        W_set_temp = [];
        fprintf(' %d.',i);
        subdata = Data_matrix(:,i-history:i-1);
        subsample_cov = cov(subdata');
        subdata = subdata';       
       if counth == 0
            % select for number of factors
            [f_number] = factor_number_selection(subdata',10);
            % get residuals from the first step time-varying PCA
            [~,~,Residuals] =Time_COV(subdata',0,f_number,tau,h);
            % cross-validation for sparse parameter lambda
            [PCV] = CV_for_spcov(Residuals,3,tau);
            % cross-validation for bandwidth
            [h] = CV_for_h(subdata',PCV,f_number,tau);
            h_set = [h_set;h];
            PCV_set = [PCV_set;PCV];
            f_number_set = [f_number_set;f_number];
            counth = counth +persist;
       else
            % updata the tunning parameter values annually
            if counth > 250
                counth = 0;
            else
                counth = counth+persist;
            end
       end
       try
            PDSCE_cov = Time_COV(subdata',PCV,f_number,tau,h);
       catch
            try_time =  try_time + 1;
       end
        
        subdata = subdata';
        w_opt1 = Markowitz_MVP(PDSCE_cov);%Markowitz_MVP(subsample_cov);%
        eq_w = (1/size(w_opt1,1))./ones(size(w_opt1,1),1);
        
        % to aviod look-ahead bias by dropping the missing data in advance
        ind_nan_set = [];
        for z = i :i+persist-1
            ind_nan = find(isnan(Data_matrix(:,z)));
            if ~isempty(ind_nan)
                ind_nan_set = [ind_nan_set;ind_nan];
            end
            w_opt1(ind_nan) = 0;
            eq_w(ind_nan) = 0;
            Holdingreturn = Data_matrix(:,z);
            Holdingreturn(ind_nan) = 0;
            PD_return = [PD_return;w_opt1'*Holdingreturn];
            eq_return = [eq_return;eq_w'*Holdingreturn];
        end
        W_set_temp = [W_set_temp,w_opt1,eq_w];           
        W_opt_set{W_count} = W_set_temp;
        W_count = W_count +1;
        ind_nan_set = unique(ind_nan_set);
        Data_matrix(ind_nan_set,:)=[];
        ind_nan_set = [];  %
    end
    fprintf('\n');
    S_PD = mean(PD_return)/std(PD_return);
    S_eq = mean(eq_return)/std(eq_return);
    info = [S_PD,S_eq, mean(PD_return),mean(eq_return),std(PD_return),std(eq_return)];
    toc;
end
