remove(list = ls())

###
### Parameters of the true distribution function
###
mean <- 50
std <- 20

sdlog     <- sqrt(log((mean/std)^2) + 1) 
meanlog   <- log(mean) - sdlog^2/2

### 
### Parameters of program 5 
###

a <- qlnorm(0.7,meanlog,sdlog)
eta <- dlnorm(a,meanlog,sdlog)
beta <- 1 - plnorm(a,meanlog,sdlog)
nu <- -DistributionPty::Dlnorm(a,2,meanlog,sdlog) # As second derivative of the distribution function

mu <- eta/nu
sigma <- 2*beta/nu

###
### Parameters of the news vendor problem
###
price <- 7 
cost  <- 1 
q <- seq(0, qlnorm(0.95, meanlog,sdlog), length = 50)

###
### Solve Equation (17)
### 

runFunc <- function(q){
  # For a given q, ...
  
  # Compute the first term of Equation (17)
  term11 <- DistributionPty::partialExpectationlnorm( min(a,q), meanlog,sdlog) 
  term12 <- diff(q*plnorm(c(q,a),meanlog,sdlog))*(a >= q)
  term1 <- price*(term11 + term12) - cost*q
  
  # Compute the optimal lower bound of the second term of Equation (17)
  H <- function(x) {
    M <- pmax(q-a,0)
    mx <- pmin(x,q-a)
    
    term1 <- q/2*(x-M)^2*(x >= M)
    term2 <- (1/6*(mx+a)^3 - 1/6*a^3  -x*a^2/2 + q^2/2*(x - M)*(x >= M))*(mx >= 0)
    
    output <-  price*(term1 + term2) 
    return(output)
  }
  boundTerm2 <- RobustTail::computeBound(H,mu,sigma,lambda = q/2, nu, direction = "min")$bound
 
  # Compute the overall value of Equation (17)
  bound <- term1 + boundTerm2 
  
  output <- data.frame(q = q, term1 = term1, boundTerm2 = boundTerm2, bound = bound)
  return(output)
}

library(plyr)
library(parallel)
optimBound <- ldply(mclapply(X = q, FUN = runFunc, mc.cores = 1)) # <- can increase the number of cores
save(optimBound, file =  "data/runNewsVendorLogNormal.RData")

# Solve NewVendor Problem
solveNewsVendor <- optim(fn = function(q) -runFunc(q)$bound,
                         par = 20, 
                         lower = 20, 
                         upper = 100,
                         method = "Brent") # we multiply by minus to maximize

# Compare to the analytical solution
q <- qlnorm((price - cost)/price, meanlog,sdlog)
data.frame(Algorithm = solveNewsVendor$par, Analytical = q)
