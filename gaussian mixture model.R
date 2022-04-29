library(tidyverse)



#### function FMM ####
FMM = function(data,weight_s,mu_s,sigma_s,epsilon,maxtimes_iterations){
  k = length(mu_s)
  n = length(data)
  if(is.null(weight_s)){weight_s = rep(1/k,k)}
  if(is.null(epsilon)){epsilon = 1e-6}
  if(is.null(maxtimes_iterations)){maxtimes_iterations = 100000}
  
  for(i in 1:maxtimes_iterations){ 
    sigma_s_old = sigma_s
    
    #E-step
    #w_ik
    w_ik = matrix(0, ncol = k, nrow = n) #k=2, obs=n=100
    for(j in 1:k){
      w_ik[, j] = weight_s[j] * dnorm(x, mean = mu_s[j], sd = sqrt(sigma_s[j]))
    }
    w_ik = w_ik / rowSums(w_ik)
    
    #M-step
    #N_k
    N_k = colSums(w_ik)
    weight_s = N_k/n
    
    for(j in 1:k){
      sigma_s[j] = sum( w_ik[, j] * (x - mu_s[j])^2) / N_k[j]
      mu_s[j] = sum(x * w_ik[, j]) / N_k[j]
    }
    
    if (max(abs(sigma_s_old - sigma_s)) < epsilon)
    {
      break
    }
    
  }
  
  list(mu=mu_s, sd=sqrt(sigma_s), weight=weight_s)          
}



library(lattice) #plot density package
#### Simulation data function #### 
data_ass = function(seed,obs_n,mean,sd,weight){
  set.seed(seed)
  n=obs_n
  mean_s = mean
  sd_s = sd
  k=length(mean_s)
  y = sample(c(seq(k)),size = n,replace=T,prob = weight)
  x=c()
  for (i in i:k) {
    x[which(y==i)] = rnorm(n=length(which(y==i)), mean = mean_s[i], sd=sd_s[i])
  }
  list(obs_number=n,mean=mean_s,sd=sd_s,data=x,index=y)
}



#### data example ####
A = data_ass(0413,2000,mean=c(1,3,5),sd=c(5,2,3),weight = c(0.5,0.2,0.3))
densityplot(~A$data, par.settings = list(plot.symbol = list(col = factor(A$index))))
FMM(A$data,weight_s=NULL,mu_s=c(3,6),sigma_s=c(10,1),epsilon = NULL, maxtimes_iterations = NULL)
FMM(A$data,weight_s=NULL,mu_s=c(1,2,3),sigma_s=c(1,1,1),epsilon = NULL, maxtimes_iterations = NULL)
#
B = data_ass(0414,1500,mean=c(1,7),sd=c(5,1),weight = c(0.25,0.75))
B = data_ass(0414,1500,mean=c(1,7),sd=c(1,5),weight = c(0.25,0.75))
densityplot(~B$data, par.settings = list(plot.symbol = list(col = factor(B$index))))
FMM(x,weight_s=NULL,mu_s=c(3,6),sigma_s=c(10,1),epsilon = NULL, maxtimes_iterations = NULL)

