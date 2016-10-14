remove(list = ls())
library(parallel)
library(plyr)
library(dplyr)
library(DistributionPty)
### 
### Parameters of the Pareto distribution
###
shape <- 1
scale <- 1

###
### Compute the optimal upper bound for various threshold a, and  intervals [c,d]
###
runFunc <- function(a, c, d)
{
  
  
  # Define H
  H  <- function(x){
    output <- 1/2*(x - max(c-a,0) )^2*( max(c-a,0) <= x) -  1/2*(x - max(d-a,0))^2*( max(d-a,0) <= x)
  }
  
  # Define the parameters of the program (5)
  eta   <- DistributionPty::Dpareto(a,d = 1,shape,scale)
  alpha <- 1 -  DistributionPty::ppareto(a,shape,scale)
  nu    <- -DistributionPty::Dpareto(a,d = 2,shape,scale)
  
  mu <- eta/nu
  sigma <- 2*alpha/nu
  
  # Compute the optimal upper bound
  out <- RobustTail::computeBound(H,mu,sigma, lambda = 0,   nu, direction = "max")
  output <- data.frame(a = a, c = c, d = d, bound = out$bound)
  return(output)
}

a  <- qpareto(seq(0.7,0.85, length=30),shape,scale)
c <- qpareto(seq(0.85,0.98, by = 0.01),shape,scale)
d <- qpareto(seq(0.86,0.99, by = 0.01),shape,scale)
p  <- ppareto(d,shape,scale) -  ppareto(c,shape,scale) # True values

runParam <- expand.grid(a = a, c = c) %>%
            transform(d = qpareto(ppareto(c, shape, scale) + 0.01, shape, scale))

optimBound <- with(runParam,mcmapply(FUN = runFunc, a = a, c = c, d = d, mc.cores = 1, SIMPLIFY = FALSE)) %>% ldply
# The number of cores can be increased to run the computations in parallel
save(optimBound,file = "data/runSyntheticPareto.RData")

