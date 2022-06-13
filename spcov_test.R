# This code get sparse residual covariance matrix estimator
# The only input is a ".mat" document generated from matlab, see more in Timevarying_factor_model.m
library(spcov)
library(PDSCE)
library(R.matlab)
loc_path2 = 'F:/matlab/matlab_true/bin/TV-MVP/'
Data <- readMat(paste(loc_path2,'Mat_R_temp.mat',sep = ""))
S <- cov(Data$subdata)
tau1 <- Data$tau[1]
output=pdsoft(s=S, lam= 0 ,tau = tau1, init = "soft")
Snew <- output$sigma
lam <- Data$P
tol <- 1e-6
fit <- try(mm <- spcov(Sigma=Snew, S=Snew, lambda=lam,
                       step.size=0.001,  trace=1,tol.outer=tol,n.inner.steps = 200,n.outer.steps = 200,thr.inner = 1e-3))
if("try-error" %in% class(fit))
{
  sigma <- Snew
}else
{
  sigma <- mm$Sigma
}
write.table(sigma,file = paste(loc_path2,'test.csv',sep = ""),col.names=F,row.names=F,sep=',') 


