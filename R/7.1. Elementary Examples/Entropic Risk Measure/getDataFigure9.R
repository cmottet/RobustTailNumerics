remove(list = ls())

### 
### Parameters of program 5 
###
a <- qexp(0.7) 
eta <- dexp(a) 
nu <- dexp(a) 
beta <- 1-pexp(a)

mu <- eta/nu
sigma <- 2*beta/nu

###
### Compute the optimal upper bound for different values of theta
### 
theta <- seq(0.04,2,length = 50)

runFunc <- function(theta){
  # Compute the first term in Equation (15)
  term1 <- 1/(theta+1)*(1-exp(-(theta+1)*a))
  
  # Compute the second term in Equation (15)
  H <- function(x) output <- exp(-theta*a)/theta^2*(theta*x + exp(-theta*x) -1)
  boundTerm2 <- RobustTail::computeBound(H,mu,sigma,lambda = 0, nu)$bound
  
  # Compute the value of Equation (15)
  bound <- 1/theta*log(term1 + boundTerm2)
  
  output <- data.frame(theta = theta, term1 = term1, boundTerm2 = boundTerm2, bound = bound)
  return(output)
}

library(plyr)
library(parallel)
optimBound <- ldply(mclapply(X = theta, FUN = runFunc, mc.cores = 1)) # Can change the number of cores! to run in parallel
save(optimBound,file = "data/runEntropyExpDist.RData")



