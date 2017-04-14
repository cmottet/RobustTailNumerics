remove(list = ls())
library(parallel)
library(plyr)
library(dplyr)
library(DistributionPty)
library(RobustTail)
### 
### Parameters of the gamma distribution
###
shape <- 2
rate <- 1

###
### Compute the optimal upper bound for various threshold a, and  intervals [c,d]
###
runFunc <- function(a, c, d)
{
  
  
  # Define H
  H  <- function(x){
    output <- 1/2*(x - max(c-a,0) )^2*( max(c-a,0) < x) -  1/2*(x - max(d-a,0))^2*( max(d-a,0) <= x)
  }
  
  # Define the parameters of the program (5)
  alpha <- 1 -  DistributionPty::Dgamma(a,D = 0, shape,rate)
  eta   <- DistributionPty::Dgamma(a,D = 1,shape,rate)
  nu    <- -DistributionPty::Dgamma(a,D = 2,shape,rate)
  
  mu <- eta/nu
  sigma <- 2*alpha/nu
  
  # Compute the optimal upper bound
  out <- RobustTail::computeBound(H,mu,sigma, lambda = 0,   nu, direction = "max")
  output <- data.frame(a = a, c = c, d = d, bound = out$bound)
  return(output)
}

a  <- qgamma(seq(0.7,0.85, length=30),shape,rate)
c <- qgamma(seq(0.85,0.98, by = 0.01),shape,rate)
d <- qgamma(seq(0.86,0.99, by = 0.01),shape,rate)
p  <- pgamma(d,shape,rate) -  pgamma(c,shape,rate) # True values

runParam <- expand.grid(a = a, c = c) %>%
            transform(d = qgamma(pgamma(c, shape, rate) + 0.01, shape, rate))

optimBound <- with(runParam,mcmapply(FUN = runFunc, a = a, c = c, d = d, mc.cores = 2, SIMPLIFY = FALSE)) %>% ldply
# The number of cores can be increased to run the computations in parallel
#save(optimBound,file = "data/runSyntheticgamma.RData")
save(optimBound,file = "R/7.1. Elementary Examples/Tail Interval Probability/Gamma/runSyntheticGamma.RData")

