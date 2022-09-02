library(MASS)
library(mvtnorm)

#moni
set.seed(220901)
mean1 = c(1,2,3)
mean2 = c(0,0,2)
sd1 = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3,ncol=3)
sd2 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
moni1 = mvrnorm(100,mean1,sd1)
moni2 = mvrnorm(100,mean2,sd2)
moni=rbind(moni1,moni2)
epsilon=0.0001
#
dat=moni

#
k=2
alpha=c(0.7,0.3)
w_ik0 = matrix(0, ncol = k, nrow = nrow(dat)) 
w_ik_old=matrix(rep(c(1,0),nrow(dat)),nrow(dat),k,byrow = T)
maxtimes_iterations=10000
mu_s1=rep(1, ncol(dat)) #2 initial mean must be different
mu_s2=rep(0, ncol(dat))
mu_s=list(mu_s1,mu_s2)
sigma_s1=diag(ncol(dat))
sigma_s2=diag(ncol(dat))
sigma_s=list(sigma_s1,sigma_s2)

#EM
GMM = function(dat,w_ik0,alpha,mu_s,sigma_s,epsilon,maxtimes_iterations){
  
for(i in 1:maxtimes_iterations){ 
  
  #E-step
  #w_ik
  for(j in 1:k){
    w_ik0[,j] = alpha[j] * dmvnorm(dat,mu_s[[j]],sigma_s[[j]])
  }
  w_ik = w_ik0 / rowSums(w_ik0)
  
  #M-step
  #N_k
  N_k = colSums(w_ik)
  alpha = N_k/nrow(dat)
  
  for(j in 1:k){
    mu_s[[j]] = (t(as.matrix(w_ik[, j]))%*%dat) / N_k[j]
    sigma_s[[j]] = alpha[j]*t((dat - matrix(rep(mu_s[[j]],nrow(dat)),nrow(dat),3,byrow = T))) %*% (dat - matrix(rep(mu_s[[j]],nrow(dat)),nrow(dat),3,byrow = T)) / N_k[j]
  }
  
  if ( sum(abs(w_ik_old - w_ik)) < epsilon)
  {
    break
  }
  
  w_ik_old=w_ik
  
}
  
  list(mu=mu_s, COV=sigma_s, weight=w_ik) 
}
  

  
ye = GMM(moni,w_ik0,alpha,mu_s,sigma_s,epsilon,maxtimes_iterations)
  
  
  








