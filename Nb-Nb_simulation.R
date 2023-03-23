
# Simulation: R code for article:
# B.Kowalczyk, R.Wieczorkowski, New Improved Poisson And Negative Binomial Item Count Techniques For Eliciting Truthful  Answers To Sensitive Questions,
# *STATISTICS IN TRANSITION new series*, March 2022, Vol. 23 No. 1, pp. 75–88, DOI 10.21307/stattrans-2022-005

# model: negative binomial - negative binomial
source("Nb-Nb_code.R")

source("Nb_code_Tian.R") # original Tian model



symulNbNb<-function(Niter,n,pi_true,p_true,r_true)
# simulation for given parameters: n,pi, r, p
{
 
  #n=250;pi_true=0.1; p_true=0.4; r_true=2
  cat("n=",n,"pi=",pi_true,"p = ",p_true,"r=",r_true,"\n")
  tab<-matrix(0,Niter,4)
  
  for (i in 1:Niter)
  {
    #cat('Iteration ',i,'\n')
    
    ## Generate sample data
    n1<-n
    n2<-n
    
    x1<-rnbinom(n1+n2,size=r_true,prob=1-p_true)
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
    
    r1_mm<-(mean(x1)^2)/(var(x1)-mean(x1))
    #info<-'0'
    if (r1_mm<0) 
    {
      cat("Parametr r <0 !!!","\n"); 
      r1_mm<-1.0001 ;
      #info<-'1';
      #i<-i-1
    }
    p1_mm<-mean(x1)/(r1_mm+mean(x1))  
                                                   
    r2_mm<-(mean(x2)^2)/(var(x2)-mean(x2))
    #info<-'0'
    if (r2_mm<0) 
    {
      cat("Parametr r <0 !!!","\n"); 
      r2_mm<-1.0001 ;
      #info<-'1';
      #i<-i-1
    }
    p2_mm<-mean(x2)/(r2_mm+mean(x2))    
    
    ## Set parameter estimates
    pi_init = 0.5
    #pi_init<-tab[i,1]
    
    pi_init = tab[i,1]  # MM
    pi_init = ifelse(pi_init<=0 | pi_init>=1,0.5,pi_init)
    
    p1_init = p1_mm
    p2_init = p2_mm
    
    ## Run EM Algorithm
    est_em <- em.algo.nbnb(x1,x2,y1,y2,pi_init,p1_init,p2_init,r1_mm,r2_mm)
    #print(est_em)
    tab[i,2]<-est_em[[1]] #ML
    
    
    ## for original method from Tian
    tab[i,3] =  mean(y2)-mean(x1) # MM0
    pi_init = tab[i,3]  
    pi_init = ifelse(pi_init<=0 | pi_init>=1,0.5,pi_init)
    
    ## Run EM Algorithm
    est_em <- em.algo.nb(x1,x2,y2,pi_init,p1_init,r1_mm)
    #print(est_em)
    tab[i,4]<-est_em[[1]] # ML0
    
    
  }
  
  return(data.frame(n=rep(n,Niter),pi=rep(pi_true,Niter),
                    p=rep(p_true,Niter), r=rep(r_true,Niter),
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
pi_list <- c(5)/100
p_list <- c(0.4,0.5,0.6)
r_list <- c(2,3)

#params<-as.matrix(expand.grid(n=n_list,p=pi_list,lambda1=lambda_list,lambda2=lambda_list))
#params<-as.matrix(expand.grid(n=n_list[[3]],p=pi_list[[4]],lambda=lambda_list[[1]],a=a_list))
#params<-params[params[,3]<=params[,4],]

params<-as.matrix(expand.grid(n=n_list,pi=pi_list,p=p_list,r=r_list))


tab <- foreach(i=1:nrow(params), .combine=rbind) %dopar% 
  {
    symulNbNb(Niter,params[i,1],params[i,2],params[i,3],params[i,4]) 
  }


t2<-Sys.time()
cat("Elapsed time     ",(t2-t1),"\n")
stopCluster(cl)


saveRDS(tab,"Nb-Nb_symul_a.rds")



