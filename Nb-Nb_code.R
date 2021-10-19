
# Maximum likelihood estimation via EM algorithm
# model: negative binomial - negative binomial

### section for X - negative binomial distribution

estep.nbnb <- function(y1,y2,pp,p1,p2,r1,r2)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group 
#   y2 - sample from second treatment group 
#   pp - estimate for parameter pi
#   p - estimate for parameter p
#   r - estimate for parameter r
#
# Values:
#   output data frame contains two columns with expected values 
#

{
  z1_estep <- (pp*y1)/(pp*y1 + (y1+r2-1)*p2*(1-pp) )
  z2_estep <- (pp*y2)/(pp*y2 + (y2+r1-1)*p1*(1-pp) )
  
  return(list(z1_estep,z2_estep))
}


mstep.nbnb <- function(x1,x2,y1,y2,e.step,r1,r2)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group 
#   y2 - sample from second treatment group
#   e.step - values from Expectation step
#
# Values:
#   output list contains values of parameters (pi,lambda,p) 
#   from M-step of EM algorithm
{
  
  # estimate pi
  pi_temp <- mean(unlist(e.step))
  
  p1_temp <-(sum(x1)+sum(y2-e.step[[2]]))/
    ( r1*(length(y1)+length(y2))+ sum(x1)+sum(y2-e.step[[2]])  ) 
  
  p2_temp <-(sum(x2)+sum(y1-e.step[[1]]))/
       ( r2*(length(y1)+length(y2))+ sum(x2)+sum(y1-e.step[[1]])  ) 
  
  
  
  
  list(pi_temp,p1_temp,p2_temp)   
}



em.algo.nbnb <- function(x1,x2,y1,y2,pi_init=NULL,p1_init=NULL,p2_init=NULL,r1=1,r2=1,
                        maxit=10000,tol=1e-6,info=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   y1 - sample from first treatment group 
#   y2 - sample from second treatment group 
#   pi_init - initial value for probability pi
#   p_init - initial vlues for parameter p
#   r - parameter r
#   maxit - maximal number of iterations (default 10000)
#   tol - values for testing convergence of algorithm (default 1-6)
#   info - logical value for printing information about iterations
#
# Values:
#   output list contains ML estimators of parameters (pi,lambda) 
# 
{
  # Initial parameter estimates
  flag <- 0
  
  if (is.null(pi_init))
  {
    if (info) cat("Initialization of parameters by method of moments ...\n")
    pi_init<-(mean(y2)-mean(x1))
    pi_init <- ifelse(pi_init<0 | pi_init>1,0.5,pi_init)
    lambda_init<-mean(x1)
    p1_init<-mean(x1)/(r1+mean(x1))  
    p2_init<-mean(x2)/(r2+mean(x2))  
    
    if (info) cat("pi_init=",pi_init," p1_init=",p1_init," p2_init=",p2_init,"\n")
  }
  
  pi_cur <- pi_init;  p1_cur <- p1_init; p2_cur <- p2_init;
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,p1_cur,p2_cur)
    new <- mstep.nbnb(x1,x2,y1,y2,estep.nbnb(y1,y2, pi_cur, p1_cur,p2_cur,r1,r2),r1,r2)
    pi_new <- new[[1]];  p1_new <- new[[2]]; p2_new <- new[[3]];
    
    #if (is.na(pi_new)) pi_new<-pi_cur
    #if (is.na(lambda_new)) lambda_new<-lambda_cur
    
    new_step <- c(pi_new,p1_new,p2_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new;  p1_cur <- p1_new; p2_cur <- p2_new;
  }
  if(!flag) warning("Did not converge\n")
  if (info) cat(i," iterations \n")
  
  return(list(pi=pi_cur, p1=p1_cur, p2=p2_cur ))
}



