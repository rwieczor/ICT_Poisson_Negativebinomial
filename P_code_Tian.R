
# Maximum likelihood estimation via EM algorithm
# based on Tian et al. (2017) paper
# model: Poisson 



estep.p <- function(y2,pp,lambda)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y2 - sample Y from second treatment group 
#   pp - estimate for parameter pi
#   lambda - estimate for parameter lambda
#
# Values:
#   output data frame contains column with expected values 
#

{
  z2_estep <- (pp*y2)/(pp*(y2)+(lambda)*(1-pp))
  
  return(list(z2_estep))
}


mstep.p <- function(x1,x2,y2,e.step)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   x1 - sample X from first treatment group
#   x2 - sample X from second treatment group
#   y2 - sample from second treatment group
#   e.step - values from Expectation step
#
# Values:
#   output list contains values of parameters (pi,lambda) 
#   from M-step of EM algorithm
{
  
  # estimate pi
  pi_temp <- mean(unlist(e.step))
  
  # estimate lambda
  lambda_temp <-(sum(x1)+sum(y2-e.step[[1]]))/(length(x1)+length(y2)) 

  list(pi_temp,lambda_temp)   
}



em.algo.p <- function(x1,x2,y2,pi_inits=NULL,lambda_inits=NULL,
                        maxit=10000,tol=1e-6,info=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   x1 - sample X from first treatment group
#   x2 - sample X from second treatment group
#   y2 - sample Y from second treatment group 
#   pi_inits - initial value for probability pi
#   lambda_inits - initial vlues for parameter lambda1, lambda2
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
  
  if (is.null(pi_inits))
  {
    if (info) cat("Initialization of parameters by method of moments ...\n")
    pi_inits<-(mean(y2)-mean(x1))
    pi_inits <- ifelse(pi_inits<0 | pi_inits>1,0.5,pi_inits)
    lambda_inits<-mean(x1)
    
    if (info) cat("pi_inits=",pi_inits," lambda_inits=",lambda_inits,"\n")
  }
  
  pi_cur <- pi_inits; lambda_cur <- lambda_inits; 
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,lambda_cur)
    new <- mstep.p(x1,x2,y2,estep.p(y2, pi_cur, lambda_cur))
    pi_new <- new[[1]]; lambda_new <- new[[2]];
    #if (is.na(pi_new)) pi_new<-pi_cur
    #if (is.na(lambda_new)) lambda_new<-lambda_cur
    
    new_step <- c(pi_new,lambda_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new; lambda_cur <- lambda_new;
  }
  if(!flag) warning("Did not converge\n")
  if (info) cat(i," iterations \n")
  
  return(list(pi=pi_cur, lambda=lambda_cur))
}



