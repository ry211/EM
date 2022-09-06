library(MASS)
library(mvtnorm)

############### moni—parameter / true value: meann,sdd################
mean1 = c(1,2,3)
mean2 = c(0,0,2)
meann=list(mean1,mean2)
sd1 = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3,ncol=3)
sd2 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sdd = list(sd1,sd2)
epsilon=0.0001
k=2
maxtimes_iterations=10000
alp = list()



################### EM-function #####################
GMM = function(dat,w_ik0,alpha,mu_s,sigma_s,epsilon,maxtimes_iterations){
  
  for(i in 1:maxtimes_iterations){ 
    
    k=ncol(w_ik0)
    
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
    alp[[i]] = alpha
    
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
  
  list(mu=mu_s, COV=sigma_s, weight=w_ik,alpha=alp,ite_time=i) 
}



############### simulation1 ######################
ye1=list()
#
set.seed(220901)
#
for (ite in 1:100) {
  
  moni1 = mvrnorm(100,mean1,sd1)
  moni2 = mvrnorm(100,mean2,sd2)
  dat  = rbind(moni1,moni2)  #means true value of alpha is c(0.5,0.5)
  alpha=c(0.7,0.3)
  w_ik0 = matrix(0, ncol = k, nrow = nrow(dat)) 
  w_ik_old=w_ik0
  
  # initial value
  mu_s1=rep(1, ncol(dat)) #2 initial mean must be different
  mu_s2=rep(0, ncol(dat))
  mu_s=list(mu_s1,mu_s2)
  sigma_s1=diag(ncol(dat))
  sigma_s2=diag(ncol(dat))
  sigma_s=list(sigma_s1,sigma_s2)
  
  #
  ye1[[ite]] = GMM(dat,w_ik0,alpha,mu_s,sigma_s,epsilon,maxtimes_iterations)
}


################# evaluation1 ####################
ye_mean1 = list()
ye_alpha1 = list()

for (ite in 1:100) {
  ye_mean1[[ite]]=matrix(unlist(ye1[[ite]]$mu),2,3,byrow = T) - matrix(unlist(meann),2,3,byrow = T)
  ye_alpha1[[ite]]=matrix(unlist(ye1[[ite]]$alpha),length(unlist(ye1[[ite]]$alpha))/2,2,byrow = T) - matrix(c(0.5,0.5),length(unlist(ye1[[ite]]$alpha))/2,2,byrow = T)
}  

avg_bias = sum(abs(unlist(ye_mean1)))/length(unlist(ye_mean1)) #0.265304
avg_bias1 = sum(unlist(ye_mean1))/length(unlist(ye_mean1)) #-0.1221698
MSE=sum((unlist(ye_mean1))^2)/length(unlist(ye_mean1)) #0.1183009
avg_alpha = sum(abs(unlist(ye_alpha1)))/length(unlist(ye_alpha1)) #0.159538








