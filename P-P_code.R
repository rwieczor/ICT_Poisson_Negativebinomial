
# Maximum likelihood estimation via EM algorithm
# model: Poisson - Poisson



estep.pp <- function(y1,y2,pp,lambda1,lambda2)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group 
#   y2 - sample from second treatment group 
#   pp - estimate for parameter pi
#   lambda1 - estimate for parameter lambda1
#   lambda2 - estimate for parameter lambda2
#
# Values:
#   output data frame contains two columns with expected values 
#

{
  z1_estep <- (pp*y1)/(pp*(y1)+(lambda2)*(1-pp))
  z2_estep <- (pp*y2)/(pp*(y2)+(lambda1)*(1-pp))
  
  return(list(z1_estep,z2_estep))
}


mstep.pp <- function(x1,x2,y1,y2,e.step)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group 
#   y2 - sample from second treatment group
#   e.step - values from Expectation step
#   a - parameter
#
# Values:
#   output list contains values of parameters (pi,lambda1,lambda2) 
#   from M-step of EM algorithm
{
  
  # estimate pi
  pi_temp <- mean(unlist(e.step))
  
  # estimate lambdas
  lambda1_temp <-(sum(x1)+sum(y2-e.step[[2]]))/(length(y1)+length(y2)) 
  lambda2_temp <-(sum(x2)+sum(y1-e.step[[1]]))/(length(y1)+length(y2)) 
  
  
  list(pi_temp,lambda1_temp,lambda2_temp)   
}



em.algo.pp <- function(x1,x2,y1,y2,pi_inits=NULL,lambda1_inits=NULL,lambda2_inits=NULL,
                        maxit=10000,tol=1e-6,info=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   y1 - sample from first treatment group 
#   y2 - sample from second treatment group 
#   pi_inits - initial value for probability pi
#   lambda_inits - initial vlues for parameter lambda1, lambda2
#   a - additional parameter of the method (multiplier for Z)
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
    lambda1_inits<-mean(x1)
    lambda2_inits<-mean(x2)
    
    if (info) cat("pi_inits=",pi_inits," lambda_inits=",lambda1_inits,lambda2_inits,"\n")
  }
  
  pi_cur <- pi_inits; lambda1_cur <- lambda1_inits; lambda2_cur <- lambda2_inits;
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,lambda1_cur,lambda2_cur)
    new <- mstep.pp(x1,x2,y1,y2,estep.pp(y1,y2, pi_cur, lambda1_cur,lambda2_cur))
    pi_new <- new[[1]]; lambda1_new <- new[[2]]; lambda2_new <- new[[3]];
    
    #if (is.na(pi_new)) pi_new<-pi_cur
    #if (is.na(lambda_new)) lambda_new<-lambda_cur
    
    new_step <- c(pi_new,lambda1_new,lambda2_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new; lambda1_cur <- lambda1_new; lambda2_cur <- lambda2_new;
  }
  if(!flag) warning("Did not converge\n")
  if (info) cat(i," iterations \n")
  
  return(list(pi=pi_cur, lambda1=lambda1_cur, lambda2=lambda2_cur))
}



