remove(list = ls())

library(dplyr)
library(plyr)
library(parallel)

# Compute the optimal upper bound for various values of eta, nu, beta
load ("data/syntDatalogNormalCI.RData")

# Get the bootstrapped sample obtained to build the CI's
# and compute point estimates
bootSample <- CI[[55]]$bootSample
Fhat <- apply(bootSample, 2,mean)

etaCI <- CI[[55]]$hyperrectangle$d1
betaCI <- CI[[55]]$hyperrectangle$m0
nuCI <- -CI[[55]]$hyperrectangle$d2

parameters <- expand.grid(eta = c(etaCI, Fhat[2]), 
                          beta = c(betaCI, Fhat[3]), 
                          nu = c(rev(nuCI), -Fhat[1])) %>%
  transform(mu = eta/nu,  sigma =  2*beta/nu) %>%
  transform(feasible = mu^2 <= sigma)

# Parameters of the problem
c <- 4
d <- 5
a <- CI[[55]]$a# A index is 55 (compare to the full vector a)  # See Analysis Synthetic Log Normal to get this
H <- function(x) output <- 1/2*( (x+a -c)^2*(x+a>=c) - (x+a -d)^2*(x+a>=d))*(x >=a)

runFunc <- function(nu,mu,sigma){
  bound <- RobustTail::computeBound(H, mu, sigma, lambda = 0, nu = nu)$bound
  output <- data.frame(bound = bound)
  return(output)
}

optimBound <- ldply(with(parameters, mcmapply(FUN = runFunc, 
                                              mu = mu,
                                              sigma = sigma,
                                              nu = nu, 
                                              mc.cores = 1, # The number of cores can be increases to run computation in parallel
                                              SIMPLIFY = FALSE)))
# Create Table 2 in Latex exportable format
library(xtable)
typeValue <- c("Lower bound", "Upper bound", "Estimated value")
dataTable <-  expand.grid(eta = typeValue, 
                          beta = typeValue, 
                          nu =   typeValue) %>%
  transform(bound = optimBound) %>%
  filter(is.finite(bound)) %>%
  arrange(bound)

names(dataTable) <-c("$ \\eta$", "$ \\beta$", "$ \\nu$", "Optimal upper bound")
caption <- "Sensitivity analysis of the optimal upper bound of $P(4<X<5)$ for the synthetic data in Example \\ref{example:synthetic}."
label <- "Tab:Ex SyntD  sensitivity"
display <- c("d","s","s","s","E")
align  <- rep("c",length(display))

xtable(dataTable, caption = caption , label = label, display = display, align = align) %>% 
  print(include.rownames=FALSE, type = "latex", sanitize.text.function = identity, file = "tables/logNormalSensitivity.tex")

