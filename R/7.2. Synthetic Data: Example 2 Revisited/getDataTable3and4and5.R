remove(list=ls())

###
### Simulate Data
###
n <- 200
N <- 100
meanlog <- 0   
sdlog <- 0.5

set.seed(50)
sample   <- matrix(rlnorm(n*N,meanlog,sdlog), ncol = n, nrow = N)
save(sample, file = "data/syntDatalogNormalProbCov.RData")

### 
### Compute the 95% confidence intervals of eta, nu, and beta
### for all 100 simulated sample, when a = 2.8 and a = 3.1
###
a <- c(2.8,3.1)

CI <- vector("list",length(a))
for (i in 1:length(a))
  CI[[i]] <- apply(sample,1,function(oneSample)
    RobustTail::getCIMomentAndDerivatives(oneSample,
                                   a[i],
                                   m = 0,
                                   d= 1:2,
                                   nboot = 1E3,
                                   alpha = 0.05,
                                   bootSample = TRUE,
                                   mc.cores = 1) # Can be increased to run computations in parallel
  )

save(CI,file = "data/syntDatalogNormalCIProbCov.RData")

###
### Get data Table 3 and 4
###
c <- seq(4,9)
d <- c + 1

optimBound <- NULL
runCI <- c(CI[[1]], CI[[2]])
for (i in 1:length(runCI)){
  
  a <- runCI[[i]]$a
  estimations <- runCI[[i]]$hyperrectangle
  
  nu <- -estimations$d2[1]
  eta <- estimations$d1
  beta <-  estimations$m0
  
  mu <- eta/nu
  sigma <- 2*beta/nu
  
  
  runFunc <-function(c, d)
  {
    
    mL <- max(c - a, 0)
    mU <- max(d - a, 0)
    
    H <- function(x) output <- 1/2*((x - mL)^2*(x>= mL) - (x - mU)^2*(x>= mU))*(x >= a)
    bound <- RobustTail::computeBound(H,mu,sigma,lambda = 0,nu)$bound
    
    sampleIndex <- if (i %% N == 0) i - (i%/%N-1)*N else i %% N
    output <- data.frame(sampleIndex  = sampleIndex,a = a, c = d, d = d, bound = bound)
    return(output)
  }
  # mc.cores can be increased to run computations in parallel
  newBounds <- parallel::mcmapply(FUN = runFunc, c = c, d = d, mc.cores = 1, SIMPLIFY = FALSE)
  optimBound <- rbind(optimBound, plyr::ldply(newBounds))
}

save(optimBound, file = "data/runSyntheticLogNormalProbCov.RData")


###
###  Get data Table 5
###
library(DistributionPty)
library(dplyr)
GPDBound <- expand.grid(sampleIndex = 1:N, c = seq(4,9)) %>% transform( d = c + 1)

runFuncGPD <- function(data, sample){
  i <- data$sampleIndex
  d <- data$d
  c <- data$c
  u <- 1.8
  
  fitGPD <- tryCatch(
    QRM::fit.GPD(sample[i,],threshold = u,type = "ml"), 
    error = function(e)e)
  
  bound <- if ("error" %in% class(fitGPD)) NA else{
    h <- function(xi, beta) QRM::pGPD(d -u,xi,beta) - QRM::pGPD(c -u,xi,beta)
    hGrad <- function(xi,beta) gradientpGPD(d-u,xi,beta) - gradientpGPD(c-u,xi,beta)
    asymptoticCIforGPDfit(fitGPD,h,hGrad,verbose = FALSE)$uB
  }
  
  output <- data.frame(bound = bound)
  return(output)
}

library(plyr)
GPDBound <- ddply(GPDBound, .(sampleIndex, c), .fun = runFuncGPD, sample = sample)
save(GPDBound, file = "data/runSyntheticLogNormalProbCovGPD.RData")

