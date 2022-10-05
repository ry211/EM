library(tidyverse)
library(MASS)
library(mvtnorm)
library(stats)
library(grid)
library(graphics)
library(ggthemes)
library(gtable)
library(gridExtra)
library(ggplot2)
library(cowplot)
seed = 2022

#### data cleaning ####
dat0= read.delim("/Users/ry/Desktop/work/practicum/application/data_practicum.txt", header=TRUE)
colnames(dat0)
dat_est = transmute(dat0, physical1=physical_promis_score_AS, social1=social_promis_score_AS, cogfunctioning1=cogfunctioning_promis_score_AS,
                    sleepdisturb1=sleepdisturb_promis_score_AS, anxiety1=anxiety_promis_score_AS, depression1=depression_promis_score_AS, 
                    fatigue1=fatigue_promis_score_AS, pain1=pain_promis_score_AS,
                    physical2=f_physical_promis_score_AS, social2=f_social_promis_score_AS, cogfunctioning2=f_cogfunctioning_promis_score_AS,
                    sleepdisturb2=f_sleepdisturb_promis_score_AS, anxiety2=f_anxiety_promis_score_AS, depression2=f_depression_promis_score_AS, 
                    fatigue2=f_fatigue_promis_score_AS, pain2=f_pain_promis_score_AS)
dat_est1 = na.omit(dat_est)
nrow(dat_est)-nrow(dat_est1) #0 - no missing here
dat_est_dif = transmute(dat_est1, physical=physical2-physical1, social=social2-social1, cogfunctioning=cogfunctioning2-cogfunctioning1,
                        sleepdisturb=sleepdisturb2-sleepdisturb1, anxiety=anxiety2-anxiety1, depression=depression2-depression1,
                        fatigue=fatigue2-fatigue1, pain=pain2-pain1)
colnames(dat_est_dif)
#
dat = matrix(unlist(dat_est_dif),2850,8)
dim(dat)

#box_plot
#https://www.r-bloggers.com/2016/01/hierarchical-clustering-in-r-2/
clusters = hclust(dist(dat))
plot(clusters)
clusterCut = cutree(clusters, 4)
table(clusterCut, final[,4]) #use k=4 as final class here
#
clusters = hclust(dist(dat), method = 'average')
plot(clusters)


############### initial value 3 ##################
initial_value = function(k,dat){
  alpha=rep(1/k,k)
  mu0 = list()
  sigma0 = list()
  for (i in 1:k) {
    mu0[[i]] = colMeans(dat)+rnorm(8,0,1)
    mu0[[i]][1:6] = 0
    #sigma0[[i]] = diag(ncol(dat))*20 #ini3(1)
    sigma0[[i]] = var(dat) #ini3(2)
  }
  list(alpha=alpha,mu0=mu0,sigma0=sigma0)
}


############### GMM 2 ##################
GMM = function(k,dat,epsilon,maxtimes_iterations){
  
  #initial value
  w_ik0 = matrix(1/(k*nrow(dat)), ncol = k, nrow = nrow(dat)) #1/k
  w_ik_old = w_ik0
  ini = initial_value(k,dat)
  alpha = ini$alpha
  mu_s = ini$mu0
  sigma_s = ini$sigma0
  ll2_old=0
  ll_ite = list()
  
  for(i in 1:maxtimes_iterations){ 
    
    #E-step
    #w_ik
    for(j in 1:k){
      w_ik0[,j] = alpha[j] * dmvnorm(dat,mu_s[[j]],sigma_s[[j]])
    }
    w_ik = w_ik0 / rowSums(w_ik0) # col-class, row-obs   
    
    #M-step
    #N_k
    N_k = colSums(w_ik,na.rm = T)
    alpha = N_k/nrow(dat)
    
    for(j in 1:k){
      mu_s[[j]] = (t(as.matrix(w_ik[, j])) %*% dat) / N_k[j]
      xx= dat - matrix(rep(mu_s[[j]],nrow(dat)),nrow(dat),ncol(dat),byrow = T)
      xx1=xx
      for (m in 1:nrow(dat)) {
        xx1[m,]=xx[m,]*w_ik[m,j]
      }
      sigma_s[[j]] = (t(xx1) %*% xx) / N_k[j]
    }
    
    ll0=list()
    for (j in 1:k) {
      ll0[[j]] = alpha[j]*dmvnorm(dat,mu_s[[j]],sigma_s[[j]])
    }
    ll1 = mapply(`+`, ll0)
    ll2 = sum(log(rowSums(ll1)))
    ll_ite[[i]]=ll2
    
    sum(abs(ll2_old - ll2) ) 
    
    if ( i>=15 &i<=25)
    {
      print(i)
      print(alpha)
      print(mu_s)
      print(sigma_s)
      print(ll2)
      print("_________________________________________")

    }
    
    if ( sum(abs(ll2_old-ll2)) < epsilon)
    {
      break
    }
    
    ll2_old = ll2   
  }
  
  list(mu=mu_s, COV=sigma_s, weight=w_ik,alpha=alpha,ite_time=i,ll=ll2,ll_ite=ll_ite) 
}


####  run GMM  ####

epsilon=0.1    #   0.001
maxtimes_iterations=10000

k=4
set.seed(seed)
ye = GMM(4,dat,epsilon,maxtimes_iterations)
ggplot(data.frame(value=unlist(ye$ll_ite), ite_time=1:length(unlist(ye$ll_ite))), aes(x = ite_time, y = value)) +
  geom_line()

ye1 = list()
res = matrix(0,8,3)
set.seed(seed)
#19
#20
#20200926

for (k in 2:9) {
  ye1[[k]] = GMM(k,dat,epsilon,maxtimes_iterations)
  ll =  ye1[[k]]$ll
  n_para = length(unlist(ye1[[k]]$mu)) + length(unlist(ye1[[k]]$COV)) + k-1 
  aic = -2*ll+2*n_para
  bic = -2*ll+log(nrow(dat))*n_para
  res[k-1,] = cbind(aic,bic,k=k)
}
colnames(res) = c("AIC","BIC","k")


# plot AIC BIC
re = gather(as.data.frame(res),type,`index value`,-k)
ggplot(data=re,mapping=aes(x=k,y=`index value`,group=type))+
  geom_line(aes(linetype=type,color=type))+
  geom_point(aes(color=type))+
  theme(legend.position="top")

# plot: 4567
# still show class num_obs
# mean value of 8 index among every class
# give class - obs, specific class for k=2 ,, 8

#use this code to access final class in different k
final = matrix(0,nrow(dat),9)
final[,1] = 1:nrow(dat)
colnames(final) = c("obs_num","k=2","k=3","k=4","k=5","k=6","k=7","k=8","k=9")
for (k in 2:9) {
  colnames(ye1[[k]][["weight"]]) = c(1:k)
  x = ye1[[k]][["weight"]]/matrix(rep(ye1[[k]][["alpha"]],2850),2850,k,byrow = T)
  aa = apply(x, 1, function(t) colnames(x)[which.max(t)])
  final[,k] = aa
}
table(final[,2])
table(final[,3])
table(final[,4])
table(final[,5])
table(final[,6])
table(final[,7])
table(final[,8])
table(final[,9])

#connect dat with final
colnames(dat) = c("physical","social","cogfunctioning","sleepdisturb","anxiety","depression","fatigue","pain")
datt1 = cbind(as.data.frame(final),dat)
colnames(datt1)
kset = list()
kplot = list()
for (k in 2:9) {
  set = cbind(class = datt1[,k],datt1[,10:17])
  kset[[k]] = list(set,aggregate((set[,2:9]),list(class=set[,1]),mean))
  kplot[[k]] = gather(kset[[k]][[2]], idx_nm, value, -class)
}

#plot mean
plot_fun = function(i){
  p1 = ggplot(data=kplot[[i]],mapping=aes(x=idx_nm,y=value,group=class))+
    geom_line(aes(linetype=class,color=class))+
    geom_point(aes(color=class))+
    theme(legend.position="right")
  freq_plot0 = as.data.frame(table(kset[[i]][[1]][,1]))
  colnames(freq_plot0) = c("Class","Frequence")
  freq_plot1 = t(freq_plot0)
  tbl = tableGrob(freq_plot1, rows=row.names(freq_plot1), theme= ttheme_minimal())
  tbl <- gtable_add_grob(tbl,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(tb), l = 1, r = ncol(tbl))
  tbl <- gtable_add_grob(tbl,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(tbl))
  grid.arrange(p1, tbl,
               nrow=2,
               as.table=TRUE,
               heights=c(8,1)
  )
  
}

plot_fun(2)
plot_fun(3)
plot_fun(4)
plot_fun(5)
plot_fun(6)
plot_fun(7)
plot_fun(8)
plot_fun(9)


#connect original data with final
datt = cbind(obs=1:nrow(dat),dat0)


#
library(mclust)
set.seed(seed)
LPA = Mclust(dat)
LPA$G
table(LPA$classification)
#plot
BIC = mclustBIC(dat) #lower means better
plot(BIC)
#??plot for dif k
edda_class = final[,7]
edda = MclustDA(data=dat, edda_class, modelType = "EDDA")

##### kmeans ####
set.seed(seed)
ye_kmeans = list()
for (k in 2:9) {
  ye_kmeans[[k]] = kmeans(dat, k)
}
#ite
#

#### GMM function ####
#install.packages("ClusterR")
library(ClusterR)
set.seed(seed)
ye_GMMfun = list()
for (k in 2:9) {
  ye_GMMfun[[k]] = GMM(dat, k)
}
#??ll - sum?



#### data analysis ####





