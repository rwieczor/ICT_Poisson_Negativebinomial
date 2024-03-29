
# Simulation: R code for article:
# B.Kowalczyk, R.Wieczorkowski, New Improved Poisson And Negative Binomial Item Count Techniques For Eliciting Truthful  Answers To Sensitive Questions,
# *STATISTICS IN TRANSITION new series*, March 2022, Vol. 23 No. 1, pp. 75–88, DOI 10.21307/stattrans-2022-005

# model: Poisson - Poisson
source("P-P_code.R")

source("P_code_Tian.R") # original Tian model


symulP<-function(Niter,n,pi_true,lambda1_true,lambda2_true)
# simulation for given parameters: n,p,lambda1, lambda2
{
  
  #n=250;pi_true=0.1; lambda1_true=1; lambda2_true=1
  cat("n=",n,"p=",pi_true,"lambda1=",lambda1_true,"lambda2=",lambda2_true,"\n")
  tab<-matrix(0,Niter,4)
  
  for (i in 1:Niter)
  {
    #cat('Iteration ',i,'\n')
    
    ## Generate sample data
    n1<-n
    n2<-n
    
    x1<-rpois(n1+n2,lambda=lambda1_true)
    x2<-rpois(n1+n2,lambda=lambda2_true)
    
    z<-rbinom(n1+n2,1,pi_true)
    
    y1<-x2[1:n1]+z[1:n1]
    y2<-x1[(n1+1):(n1+n2)]+z[(n1+1):(n1+n2)]
    
    x1<-x1[1:n1]
    x2<-x2[(n1+1):(n1+n2)]
    
    Ahat<-var(y2)/n2+var(x1)/n1
    Bhat<-var(y1)/n1+var(x2)/n2
    wopt<-Bhat/(Ahat+Bhat)
    
    tab[i,1]<-wopt*(mean(y2)-mean(x1)) + (1-wopt)*(mean(y1)-mean(x2))  # MM
                                                   
    
    ## Set parameter estimates
    pi_init = 0.5
    #pi_init<-tab[i,1]
    lambda1_init = mean(x1)
    lambda2_init = mean(x2)
    
    pi_init = tab[i,1]  # MM
    pi_init = ifelse(pi_init<=0 | pi_init>=1,0.5,pi_init)
    
    ## Run EM Algorithm
    
    est_em <- em.algo.pp(x1,x2,y1,y2,pi_init,lambda1_init,lambda2_init)
    #print(est_em)
    tab[i,2]<-est_em[[1]] #ML
    
    
    
    ## for original method from Tian
    tab[i,3] = mean(y2)-mean(x1)  # MM0
    pi_init = tab[i,3]
    pi_init = ifelse(pi_init<=0 | pi_init>=1,0.5,pi_init)
    
    ## Run EM Algorithm
    
    est_em <- em.algo.p(x1,x2,y2,pi_init,lambda1_init)
    #print(est_em)
    tab[i,4]<-est_em[[1]] #ML0
    
    
  }
  
  return(data.frame(n=rep(n,Niter),p=rep(pi_true,Niter),
                    lambda1=rep(lambda1_true,Niter),
                    lambda2=rep(lambda2_true,Niter),
                    mm=tab[,1],ml=tab[,2],mm0=tab[,3],ml0=tab[,4]))
  
}



 


# parallel computations

library(doParallel)
library(foreach)
library(doRNG)

cl <- makePSOCKcluster(4)
registerDoParallel(cl)

t1<-Sys.time()
Niter<-10000

n_list<-c(250,500,1000)
pi_list<-c(5,10,20,30)/100
pi_list<-c(5)/100
lambda_list<-c(1,2,3)


#params<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=lambda_list,lambda2=lambda_list))
#params<-as.matrix(expand.grid(n=n_list[[3]],p=pi_list[[4]],lambda=lambda_list[[1]],a=a_list))
#params<-params[params[,3]<=params[,4],]

params1<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=1,lambda2=1))
params2<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=2,lambda2=2))
params3<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=3,lambda2=3))
params4<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=1,lambda2=3))
params<-rbind(params1,params2,params3,params4)

                   


tab <- foreach(i=1:nrow(params), .combine=rbind) %dopar% 
  {
    symulP(Niter,params[i,1],params[i,2],params[i,3],params[i,4]) 
  }


t2<-Sys.time()
cat("Elapsed time     ",(t2-t1),"\n")
stopCluster(cl)


saveRDS(tab,"P-P_symul_a.rds")




