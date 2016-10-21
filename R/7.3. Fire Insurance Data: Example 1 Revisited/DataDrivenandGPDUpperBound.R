remove(list = ls())

###
### Compute the optimal upper bound using program (EC.19)
###
load("data/datafireinsuranceCI.RData")

a <- CI[[length(CI)]]$a
estimations <- CI[[length(CI)]]$hyperrectangle
eta <- estimations$d1
beta <- estimations$m0

nu <- -estimations$d2[1]
mu <- eta/nu
sigma <- 2*beta/nu

c <- 50
d <- 200
H <- function(x)
{
  Mc <- max(c-a,0)
  Md <- max(d-a,0)
  
  term1 <- 1/6*((x+a-c)^3 - (Mc+a-c)^3) - 1/2*(Mc + a-c)^2*(x -Mc)
  term2 <- 1/6*((x+a-d)^3 - (Md+a-d)^3) - 1/2*(Md + a-d)^2*(x -Md)
  
  output <- term1*(x >= Mc) - term2*(x >= Md)
  return(output)
}
optimBound <- RobustTail::computeBound(H, mu, sigma, lambda = 0,nu)

###
### GPD approach
###
load("data/datafireinsurance.RData")

library(DistributionPty)
u <- 10
fitGPD <- QRM::fit.GPD(sample$Loss.in.DKM,threshold = u,type = "ml")
h <- function(xi, beta) {
  if (xi <0) return("xi < 0. This case need to be added")
  
  term1 <- (1-QRM::pGPD(d -u,xi,beta))*(d  - u + beta)/(xi - 1)
  term2 <- (1-QRM::pGPD(c -u,xi,beta))*(max(c,u) - u + beta)/(xi - 1)
  term3 <- (d - u)*(1-QRM::pGPD(d -u,xi,beta))
  term4 <- (c - u)*(1-QRM::pGPD(c -u,xi,beta))
  
  xf <- if(xi >= 0) Inf else u - beta/xi
  output <- (term1 - term2 + term3 - term4)*(min(xf,d) >= max(c,u))
  return(output)
}

hGrad <- function(xi,beta) 
{
  term11Grad <- - gradientpGPD(d -u,xi,beta)*(d + beta - u)/(xi - 1) 
  term12Grad <- matrix(c(-(d -u +beta),xi -1),nrow = 2)*(1-QRM::pGPD(d -u,xi,beta))/(xi - 1)^2
  term1Grad <- term11Grad + term12Grad 
    
  term21Grad <- - gradientpGPD(c -u,xi,beta)*(max(u,c) + beta - u)/(xi - 1) 
  term22Grad <- matrix(c(-(max(u,c) -u +beta),xi -1),nrow = 2)*(1-QRM::pGPD(max(u,c) -u,xi,beta))/(xi - 1)^2
  term2Grad <- term21Grad + term22Grad 
  
  term3Grad <- - (d - u)*gradientpGPD(d -u,xi,beta)
  term4Grad <- - (c - u)*gradientpGPD(c -u,xi,beta)
  
  xf <- if(xi >= 0) Inf else u - beta/xi
  output <- (term1Grad - term2Grad + term3Grad - term4Grad)*(min(xf,d) >= max(c,u))
  return(output)
}
GPDbound <- asymptoticCIforGPDfit(fitGPD,h,hGrad,verbose = FALSE)

###
### Compare GPD bound and optimization approach
###
data.frame(AlgoCI = optimBound$bound, GPD = GPDbound$uB)
