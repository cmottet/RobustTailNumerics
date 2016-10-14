###
### The goal of this function is to fit a mixture of Shifted Pareto and 2-PLT
### that is feasible for program (1), i.e. we find two distribution functions 
### F1 and F2 such that
### 
### w1 + w2 = 1
### w1 f1(a) + w2 f2(a) = eta
### - w1 f1'(a) - w2 f2'(a) = nu
### w1 (1 - F1(a)) + w2 (1- F2(a)) = beta
###
### where 
###  * F1(x) is a two point mass distribution function (2-PLT)
###  * F2 is a pareto distribution, i.e. F2(x) = 1 - (scale/x)^shape for all x => scale
###  and the shape is a number between 0 and 2
###
### The problem of finding such distribution functions is equivalent to finding
### F2(x) such that
### 
### eta1 = f1(a) = (eta - w2 f2(a))/w1 = (eta - w2 eta2(a))/w1  => 0
### nu1 = -f1'(a) = (nu + w2 f2'(a))/w1 = (nu - w2 nu2(a))/w1 => 0
### beta1 =  (beta - w2 (1 - F2(a)))/w1  = (beta - w2 beta2)/w1 => 0
###
### and 
###
### mu1^2 = (eta1/nu1)^2 <= sigma1 = 2*beta1/nu1 
### 
### If there exist a pareto distribution F2(x) such that the 3 previous 
### equalities and the inequality hold, then there exists a 2-PLT distribution whose
### key parameters (x1,x2) and (p1,p2), as given in Equation (6), can
### be obtained using the function "getDistribution" available in the package "RobustTail"
### 
### The process is not deterministic, it requires sampling some parameters
### from uniform distributions. In a reproducibility perspective,
### a "seed" parameter can be passed to the function to set the seed when sampling
### from the uniform distributions
###
### we now test our function 
### 
# a <- qexp(0.7)
# eta <- dexp(a)
# nu <- dexp(a)
# beta <- 1 - pexp(a)
# 
# mixture <- fitMixturePareto2PLT(a,nu,eta,beta)
# 
# # Check the system
# with(mixture, data.frame( nu = with(Pareto, w*nu) + with(PLT, w*nu),
#                           eta = with(Pareto, w*eta) + with(PLT, w*eta),
#                           beta = with(Pareto, w*beta) + with(PLT, w*beta)))
# 
# # Check that eta1, beta1, nu1, eta2, beta2, nu2 >= 0
# with(mixture, data.frame( nu = with(Pareto, nu>= 0) & with(PLT, nu >= 0),
#                           eta =with(Pareto, eta >= 0) & with(PLT, eta >= 0),
#                           beta = with(Pareto, beta>= 0) & with(PLT, beta >= 0)))
# 
fitMixturePareto2PLT <- function(a,nu,eta,beta, iterMax = 1E4,seed = NULL)
{
  mu <- eta/nu
  sigma <- 2*beta/nu
  
  for (i in 1:iterMax){
    
    # Assign randomly w2 and the shape of F2(x) 
    if (!is.null(seed)) set.seed(seed) 
    w2 <- runif(1,0, 1)
    
    if (!is.null(seed)) set.seed(seed) 
    shape <- runif(1,0,2)
    
    w1 <- (1-w2)
    
    
    # Compute the necessary (but not sufficient) bounds on the X = 
    # (x/scale)^(-shape) of F2(x) for there existence of F1(x)
    
    # XuB is an upper bound on X that ensures that 
    # eta1, nu1, and beta1 => 0
    XuB <- min(c(1,beta/w2, eta*a/(w2*shape), nu*a^2/(w2*shape*shape+1))) 
    
    # we now compute XlB, a necessary lower bound for mu1^1 <= sigma1
    b <- mu*a*shape - (shape+1)*shape*sigma/2  - a^2
    fourac <- shape*(shape+2)*a^2*(sigma-mu^2)
    delta <- b^2 - fourac
    
    XlB <- if (delta >= 0) max(c(nu/(4*w2)*(b + sqrt(delta)), 0)) else 0
    
    # Get a random value of X respecting the bounds
    # and derive the value of the scale parameter in F2(x)
    if (!is.null(seed)) set.seed(seed)
    X <- runif(1,XlB,XuB)
    scale <- X^(1/shape)*a
    
    # We now need to compute eta1, beta1, and nu1 
    # to check if there exists a 2-PLT F1(x) for
    # this specific F2(x)
    eta2 <- DistributionPty::Dpareto(a, d = 1, scale,shape)
    nu2 <- -DistributionPty::Dpareto(a, d = 2, scale,shape)
    beta2 <- 1 - DistributionPty::ppareto(a, scale,shape)
    
    eta1 <- (eta -  w2*eta2)/w1
    nu1 <- (nu -  w2*nu2)/w1
    beta1 <- (beta -  w2*beta2)/w1
    
    mu1 <- eta1/nu1
    sigma1 <- 2*beta1/nu1
    
    # Check the necesserary condition for the existence
    # of F1(x) 
    if (mu1^2 <= sigma1) break
  }
  
  if (i > iterMax) return("No solution")
  
  P <- RobustTail::getDistribution(mu1,sigma1)
  
  Pareto <- list(scale = scale, 
                shape = shape, 
                eta = eta2, 
                nu = nu2, 
                beta = beta2, 
                w = w2)
  
  PLT <- list(x = P$x,
             p = P$p,
             w = w1,  
             nu = nu1,
             eta = eta1,
             beta = beta1)
  
  output <- list(Pareto = Pareto, PLT = PLT) 
                 
  return(output)
}

