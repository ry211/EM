library(MASS) 
library(mvnfast)

#generate data
mmean0 = c(1,1,2)
mcov0 =  t(matrix(c(3,0,1, 0,2,1, 0,0,1), nrow=3, ncol=3))
da = mvrnorm(1000, mmean0, mcov0)

#function
#mean matrix
mmean = t(matrix(c(2,1,3, 1,2,3),nrow=3,ncol=2)) #k=2,n(x)=3
#Cov matrix
mcov = t(matrix(c(3,2,1, 0,2,1, 0,0,1, 1,0,0, 0,1,0, 0,0,1), nrow=3, ncol=6))

#
maxtimes_iterations=10000
k=nrow(mmean)
#

#
n=nrow(da)
weight0 = rep(1/k,k)
#

#
for(i in 1:maxtimes_iterations){ 
  w_ik_old = w_ik
  
  #E-step
  #w_ik
  w_ik = matrix(0, ncol = k, nrow = n) #k=2, obs=n=100
  for(j in 1:k){
    w_ik[, j] = weight0[j] * dmvn(da, mmean[j,] , mcov[(1:3)+(j-1)*3,]) # a_k * p_k
  }
  wik = w_ik
  wik = w_ik / rowSums(w_ik)
  
  #M-step
  #N_k
  N_k = colSums(wik)
  weight0 = N_k/n
  
  for(j in 1:k){
      mmean[j,] = colSums(wik[, j]*da) / N_k[j]
      ########################
      mcov[(1:3)+(j-1)*3,] = (wik[, j] * (t(sweep(da[,1:3],2,mmean[j,]))))%*%((sweep(da[,1:3],2,mmean[j,]))) / N_k[j]
  }

  if (abs(sum(w_ik_old) - sum(w_ik)) < epsilon) #not sure
  {
    break
  }
  
  list(  )
}

log_like = sum(log(rowSums(w_ik)))

#list(mu=mu_s, sd=sqrt(sigma_s), weight=weight_s)  
#cat(I, “th iterations, “  Log-likelihood value  “, sigma_s_old);


}


  
