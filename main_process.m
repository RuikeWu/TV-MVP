%example 1
clear;clc;
d= 50; T = 250; rho = 0.7;
h = (2.35/sqrt(12))*T^(-0.2)*d^(-0.1);

[Y,~,s_cov]=data_generate(d,T,rho);

Kernel_set = zeros(T,T);
for  r = 1 :T
    for t = 1 :T
       k = (h^-1) * Kernel2(T,t,r,h);
       Kernel_set(t,r) = k^0.5; 
    end
end
K = factor_number_selection(Y,10);
[Factor_loadings,Factor] = Timevarying_factor_model(Y,2,Kernel_set,h);

Fir_loading= [];
for i = 1 : T
   loading_temp = Factor_loadings{i};
   Fir_loading(i,:) = loading_temp(1,:);
end
plot(Fir_loading);

%%
% example 2
tau = 0.001;

[~,~,Residuals] = Time_COV(Y,0,K,tau);
lambda = CV_for_spcov(Residuals,3,tau);
[Sigma_r,Sigma_e] = Time_COV(Y,lambda,K,tau);

map = [1, 1, 0; 0, 0.8, 1; 0, 0.6, 1;
       0, 0.4, 1;0, 0.2, 1; 0, 0, 1];
subplot(1,2,1);imagesc(s_cov);colormap(map)
subplot(1,2,2);imagesc(Sigma_e);colormap(map)
colorbar
%%   
% example 4
[J_pT_set,J_pT] = FWYZtest_step1_bt(Y,K,Kernel_set);
p_val = sum(J_pT < J_pT_set)/199;
%%
% example 5
[Ln,ifa]= FWYZtest_step2(Residuals,0.05);

%%
% example 3
clear;clc;
load('example_data.mat');
[info,Tv_return,eq_return] = Rolling_window(R0,250,5,251,500);
plot(cumsum(Tv_return))
hold on ;
plot(cumsum(eq_return),'--')
legend({'Tv\_return','Equal\_return'},'Location','northwest','NumColumns',2)

%%
% example 6
alpha1 = 0.05;
d = 100;
T = 300;
loop = 500;
level = size_test(d,T,loop,alpha1);
%% Large example in section Example
clear;clc;
d = 50; T = 200; rho = 1;
h = (2.35/sqrt(12))*T^(-0.2)*d^(-0.1);

rng(2);
%generate virtual excess return data
[Y,r_cov,s_cov]=data_generate(d,T,rho); 

%with excess return data in hand, we first test for time varyingness of factor
%loading and residual covariance matrix
alpha1 = 0.05; % set the significant level
Kernel_set = zeros(T,T);
for  r = 1 :T
    for t = 1 :T
       k = (h^-1) * Kernel2(T,t,r,h);
       Kernel_set(t,r) = k^0.5; 
    end
end
K = factor_number_selection(Y,10); % estimate the number of factors by BIC
[J_pT_set,J_pT] = FWYZtest_step1_bt(Y,K,Kernel_set);
p_val = sum(J_pT < J_pT_set)/199;


[~,~,Residuals] = Time_COV(Y,0,K,0);
[Ln,ifa]= FWYZtest_step2(Residuals,alpha1);

tau = 0.001;
lambda = CV_for_spcov(Residuals,3,tau);
[Sigma,RES_cov,Residuals] = Time_COV(Y,lambda,K,tau);
sample_cov = cov(Y');

wopt  = Markowitz_MVP(r_cov);
wopt_tv = Markowitz_MVP(Sigma);
wopt_sam = Markowitz_MVP(sample_cov);

plot(wopt,'xr-'); hold on; plot(wopt_tv,'ob--');

devia_tv = abs(wopt - wopt_tv);
devia_sam = abs(wopt - wopt_sam);
plot(devia_tv,'r');hold on ; plot(devia_sam,'b--');


risk = wopt'* r_cov *wopt;
risk_tv = wopt_tv'*Sigma*wopt_tv;
risk_sam = wopt_sam'*sample_cov * wopt_sam;



