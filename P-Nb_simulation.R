
# Simulation: R code for article:
# B.Kowalczyk, R.Wieczorkowski, New Improved Poisson And Negative Binomial Item Count Techniques For Eliciting Truthful  Answers To Sensitive Questions,
# *STATISTICS IN TRANSITION new series*, March 2022, Vol. 23 No. 1, pp. 75–88, DOI 10.21307/stattrans-2022-005

# model: Poisson - negative binomial
source("P-Nb_code.R")


symulPNb<-function(Niter,n,pi_true,lambda_true,p_true,r_true)
# simulation for given parameters: n,pi,lambda, r, p
{
  
  cat("n=",n,"pi=",pi_true,"lambda=",lambda_true,"p = ",p_true,"r=",r_true,"\n")
  tab<-matrix(0,Niter,2)
  
  for (i in 1:Niter)
  {
    #cat('Iteration ',i,'\n')
    
    ## Generate sample data
    n1<-n
    n2<-n
    
    x1<-rpois(n1+n2,lambda=lambda_true)
    x2<-rnbinom(n1+n2,size=r_true,prob=1-p_true) # inna parametryzacja w R
    
    z<-rbinom(n1+n2,1,pi_true)
    
    y1<-x2[1:n1]+z[1:n1]
    y2<-x1[(n1+1):(n1+n2)]+z[(n1+1):(n1+n2)]
    
    x1<-x1[1:n1]
    x2<-x2[(n1+1):(n1+n2)]
    
    Ahat<-var(y2)/n2+var(x1)/n1
    Bhat<-var(y1)/n1+var(x2)/n2
    wopt<-Bhat/(Ahat+Bhat)
    
    tab[i,1]<-wopt*(mean(y2)-mean(x1)) + (1-wopt)*(mean(y1)-mean(x2))  # MM
                                                   
    r_mm<-(mean(x2)^2)/(var(x2)-mean(x2))
    #info<-'0'
    if (r_mm<0) 
    {
      cat("Parametr r <0 !!!","\n"); 
      r_mm<-1.00001 ;
      #info<-'1';
      #i<-i-1
    }
    p_mm<-mean(x2)/(r_mm+mean(x2))    
    
    ## Set parameter estimates
    pi_init = 0.5
    #pi_init<-tab[i,1]
    lambda_init = mean(x1)
   
    pi_init = tab[i,1]  # MM
    pi_init = ifelse(pi_init<=0 | pi_init>=1,0.5,pi_init)
    
    p_init = p_mm
    
    ## Run EM Algorithm
    est_em <- em.algo.pnb(x1,x2,y1,y2,pi_init,lambda_init,p_init,r_mm)
    #print(est_em)
    tab[i,2]<-est_em[[1]] #ML
    
    
  }
  
  return(data.frame(n=rep(n,Niter),pi=rep(pi_true,Niter),
                    lambda=rep(lambda_true,Niter),
                    p=rep(p_true,Niter), r=rep(r_true,Niter),
                    mm=tab[,1],ml=tab[,2]))
  
}






# parallel computations

library(doParallel)
library(foreach)
library(doRNG)

cl <- makePSOCKcluster(8)
registerDoParallel(cl)

t1<-Sys.time()
Niter<-10000

n_list<-c(100,250,500,1000)
pi_list<-c(5,10,15,20,25,30)/100
lambda_list<-c(1,2,3,4)


#params<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=lambda_list,lambda2=lambda_list))
#params<-as.matrix(expand.grid(n=n_list[[3]],p=pi_list[[4]],lambda=lambda_list[[1]],a=a_list))
#params<-params[params[,3]<=params[,4],]

params<-as.matrix(expand.grid(n=n_list,pi=pi_list,lambda=2,p=0.4,r=2))


tab <- foreach(i=1:nrow(params), .combine=rbind) %dopar% 
  {
    symulPNb(Niter,params[i,1],params[i,2],params[i,3],params[i,4],params[i,5]) 
  }


t2<-Sys.time()
cat("Elapsed time     ",(t2-t1),"\n")
stopCluster(cl)


saveRDS(tab,"P-Nb_symul.rds")




