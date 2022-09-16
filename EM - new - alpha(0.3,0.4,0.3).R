library(MASS)
library(mvtnorm)

############### moni—parameter / true value: meann,sdd################
mean1 = c(1,2,3)
mean2 = c(0,0,2)
mean3 = c(-2,-1,0)
meann=list(mean1,mean2,mean3)
sd1 = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3,ncol=3)
sd2 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sd3 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sdd = list(sd1,sd2,sd3)
epsilon=0.0001
k=3
maxtimes_iterations=10000


################### EM-function #####################
initial_value = function(k,dat){
  alpha=rep(1/k,k)
  mu0 = list()
  sigma0 = list()
  for (i in 1:k) {
    mu0[[i]] = rep(i, ncol(dat)) 
    sigma0[[i]] = diag(ncol(dat))
  }
  list(alpha=alpha,mu0=mu0,sigma0=sigma0)
}

#
GMM = function(k,dat,epsilon,maxtimes_iterations){
  
  #initial value
  w_ik0 = matrix(0, ncol = k, nrow = nrow(dat)) 
  w_ik_old = w_ik0
  k=ncol(w_ik0)
  ncol_dat = ncol(dat)
  ini = initial_value(k,dat)
  alpha = ini$alpha
  mu_s = ini$mu0
  sigma_s = ini$sigma0
  
  for(i in 1:maxtimes_iterations){ 
    
    #E-step
    #w_ik
    for(j in 1:k){
      w_ik0[,j] = alpha[j] * dmvnorm(dat,mu_s[[j]],sigma_s[[j]])
    }
    w_ik = w_ik0 / rowSums(w_ik0) # col-class, row-obs
    
    #M-step
    #N_k
    N_k = colSums(w_ik)
    alpha = N_k/nrow(dat)
    
    for(j in 1:k){
      mu_s[[j]] = (t(as.matrix(w_ik[, j]))%*%dat) / N_k[j]
      xx= dat - matrix(rep(mu_s[[j]],nrow(dat)),nrow(dat),3,byrow = T)
      xx1=xx
      for (m in 1:nrow(dat)) {
        xx1[m,]=xx[m,]*w_ik[m,j]
      }
      sigma_s[[j]] = t(xx1) %*% xx / N_k[j]
    }
    
    if ( sum(abs(w_ik_old - w_ik)) < epsilon)
    {
      break
    }
    
    w_ik_old=w_ik
    
  }
  
  ll0=list()
  for (j in 1:k) {
    ll0[[j]] = alpha[j]*dmvnorm(dat,mu_s[[j]],sigma_s[[j]])
  }
  ll1 = mapply(`+`, ll0)
  ll2 = sum(log(rowSums(ll1))) 
  

  list(mu=mu_s, COV=sigma_s, weight=w_ik,alpha=alpha,ite_time=i,ll=ll2) 
}



############### simulation1 ######################
ye1=list()
#
set.seed(220912)
#
for (ite in 1:100) {
  
  moni1 = mvrnorm(300,mean1,sd1)
  moni2 = mvrnorm(400,mean2,sd2)
  moni3 = mvrnorm(300,mean3,sd3)
  dat  = rbind(moni1,moni2,moni3)  #means true value of alpha is c(0.3,0.4,0.3)

  ye1[[ite]] = GMM(k,dat,epsilon,maxtimes_iterations)
}


################# evaluation1 ####################
ye_mean1 = list()
ye_alpha1 = list()
#ite=1

for (ite in 1:100) {
  ord_idx = order(rowSums(matrix(unlist(ye1[[ite]]$mu),3,3,byrow = T)),decreasing = T)
  ye_mean1[[ite]]=matrix(unlist(ye1[[ite]]$mu),3,3,byrow = T)[ord_idx,] - matrix(unlist(meann),3,3,byrow = T)
  ye_alpha1[[ite]]=sort(unlist(ye1[[ite]]$alpha))-c(0.3,0.3,0.4)}  

avg_bias = sum(abs(unlist(ye_mean1)))/length(unlist(ye_mean1)) #0.1417739
avg_bias1 = sum(unlist(ye_mean1))/length(unlist(ye_mean1)) #0.03354108
MSE=sum((unlist(ye_mean1))^2)/length(unlist(ye_mean1)) #0.05374187
avg_alpha = sum(abs(unlist(ye_alpha1)))/length(unlist(ye_alpha1)) #0.03605054
avg_alpha1 = sum(unlist(ye_alpha1))/length(unlist(ye_alpha1)) #-9.251859e-19







