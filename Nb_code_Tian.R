
# Maximum likelihood estimation via EM algorithm
# based on Tian et al. (2017) paper
# model: negative binomial


estep.nb <- function(y2,pp,p,r)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y2 - sample Y from second treatment group 
#   pp - estimate for parameter pi
#   p - estimate for parameter p
#   r - estimate for parameter r
#
# Values:
#   output data frame contains column with expected values 
#

{
  z_estep <- (pp*y2)/(pp*y2 + (y2+r-1)*p*(1-pp) )
  
  return(list(z_estep))
}


mstep.nb <- function(x1,x2,y2,e.step,r)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   x1 - sample X from first treatment group
#   x2 - sample X from second treatment group
#   y2 - sample Y from second treatment group
#   e.step - values from Expectation step
#
# Values:
#   output list contains values of parameters (pi,lambda,p) 
#   from M-step of EM algorithm
{
  
  # estimate pi
  pi_temp <- mean(unlist(e.step))
  
  p_temp <-(sum(x1)+sum(y2-e.step[[1]]))/
    ( r*(length(x1)+length(y2))+ sum(x1)+sum(y2-e.step[[1]])  ) 
  

  list(pi_temp,p_temp)   
}



em.algo.nb <- function(x1,x2,y2,pi_init=NULL,p_init=NULL,r=1,
                        maxit=10000,tol=1e-6,info=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   x1 - sample X from first treatment group
#   x2 - sample X from second treatment group
#   y2 - sample Y from second treatment group 
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
    p_init<-mean(x1)/(r+mean(x1))  
    
    if (info) cat("pi_init=",pi_init," p_init=",p_init,"\n")
  }
  
  pi_cur <- pi_init;  p_cur <- p_init;
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,p_cur)
    new <- mstep.nb(x1,x2,y2,estep.nb(y2, pi_cur, p_cur,r),r)
    pi_new <- new[[1]];  p_new <- new[[2]];
    
    #if (is.na(pi_new)) pi_new<-pi_cur
    #if (is.na(lambda_new)) lambda_new<-lambda_cur
    
    new_step <- c(pi_new,p_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new;  p_cur <- p_new;
  }
  if(!flag) warning("Did not converge\n")
  if (info) cat(i," iterations \n")
  
  return(list(pi=pi_cur, p=p_cur ))
}



