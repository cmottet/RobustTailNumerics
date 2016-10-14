remove(list = ls())

library(plyr)
library(dplyr)

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
### Create the data for the plot
###

# Get the results of Equation (15) 
load("data/runEntropyExpDist.RData")
term1 <- optimBound$term1
theta <- optimBound$theta

# Entropic risk measure value for the true standard exponential distribution
dataPlot <- select(optimBound, theta, value = bound) %>%  
  transform(method = "Optimal Upper Bound") %>%
  rbind(data.frame(theta, value = -1/theta*log(theta + 1),method = "Exponental Density"))

# Entropic risk measure value for some 2-PLT feasible solution of program (2)
P <- RobustTail::getDistribution(mu,sigma, x1 = mu/5)
term2 <- sapply(theta, function(theta) {
  H <- function(x) output <- exp(-theta*a)/theta^2*(theta*x + exp(-theta*x) -1)
  with(P, nu*(p[1]*H(x[1]) + p[2]*H(x[2])))
})

dataPlot <- dataPlot %>%
  rbind (data.frame(theta, value =  1/theta*log(term1 + term2), method = "2-PLT"))

# Entropic risk measure value for some mixture of shifted pareto and 2-PLT feasible solution of program (2)
source("R/extras.R")
mixture <- fitMixturePareto2PLT(a,nu,eta,beta, seed=200)
H <- function(x) exp(-theta*a)/theta^2*(theta*x + exp(-theta*x) -1)
term2PLT <- with(mixture$PLT, nu*(p[1]*H(x[1]) +  p[2]*H(x[2])))

library(DistributionPty)
set.seed <- 100
x <- with(mixture$Pareto, rpareto(n = 1e5, scale, shape))
grid <- expand.grid(theta = theta, x = x )
term2Pareto <- ddply(grid, .(theta), function(data) with(data, mean(exp(-theta*x)*(x >= a))))$V1
term2 <- with(mixture, PLT$w*term2PLT + Pareto$w*term2Pareto)

dataPlot <- dataPlot %>%
  rbind(data.frame(theta, value =  1/theta*log(term1 + term2), method = "Mixture of Pareto and 2-PLT"))


###
### Plot results
###
bitmap("pics/entropy_bound_exp_70th_percentile.tiff",res = 300, width = 5,height = 5) # Save in .tiff for better resolution...
library(ggplot2)
ggplot(dataPlot, aes(x = theta, y = value)) + 
  geom_line(aes(linetype = method)) + 
  labs(x = expression(theta), y = "Entropic Risk") + 
  labs(linetype = "") + 
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) 
dev.off()
# The command below only works for Mac OS X systems
# It converts to a png format without loss 
# system("sips -s format png pics/entropy_bound_exp_70th_percentile.tiff --out pics/entropy_bound_exp_70th_percentile.png")

