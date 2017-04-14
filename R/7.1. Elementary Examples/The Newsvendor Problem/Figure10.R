remove(list = ls())

library(DistributionPty)
library(plyr)
library(dplyr)
library(svglite)


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

###
### Create the data for the plot
###

# Get the results of Equation (17)  for a given q
load("data/runNewsVendorLogNormal.RData")
q <- optimBound$q
term1 <- optimBound$term1

dataPlot <- select(optimBound, q, term1, value = bound)
dataPlot <- transform(dataPlot, method = "Optimal Upper Bound")

# Equation (17)  for the true log normal distribution
term21 <- q*(1 - plnorm(pmax(q,a), meanlog,sdlog))
term22 <- (partialExpectationlnorm(q, meanlog,sdlog) -  partialExpectationlnorm(a, meanlog,sdlog))*(q >=a)
term2 <- price*(term21 + term22)
dataPlot <- rbind(dataPlot, data.frame(q, term1, value = term1 + term2, method = "Log-normal Density"))

# Equation (17)  for some 2-PLT feasible solution of program (2)
P <- RobustTail::getDistribution(mu,sigma, x1 = mu/5)
term2 <- sapply(q, function(q) {
  H <- function(x) {
    M <- pmax(q-a,0)
    mx <- pmin(x,q-a)
    
    term1 <- q/2*(x-M)^2*(x > M)
    term2 <- (1/6*(mx+a)^3 - 1/6*a^3  -x*a^2/2 + q^2/2*(x +a -q)*(x + a > q))*(mx >= 0)
    price*(term1 + term2) 
  }
  
  with(P, nu*(p[1]*H(x[1]) + p[2]*H(x[2])))
})

dataPlot <- rbind(dataPlot, data.frame(q, term1, value = term1 + term2, method = "2-PLT"))

# Equation (17) for some mixture of shifted pareto and 2-PLT feasible solution of program (2)
source("R/extras.R")
mixture <- fitMixturePareto2PLT(a,nu,eta,beta, seed=200)

#  PLT term of the objective value
H <- function(x) {
  M <- pmax(q-a,0)
  mx <- pmin(x,q-a)
  
  term1 <- q/2*(x-M)^2*(x > M)
  term2 <- (1/6*(mx+a)^3 - 1/6*a^3  -x*a^2/2 + q^2/2*(x + a -q)*(x + a >= q))*(mx >= 0)
  
  output <-  price*(term1 + term2) 
  return(output)
}
term2PLT <- with(mixture$PLT, nu*(p[1]*H(x[1]) + p[2]*H(x[2])))

# Pareto term of the objective value
set.seed <- 100 ; D <- with(mixture$Pareto, rpareto(n = 1e5, scale, shape))
grid <- expand.grid(q = q, D = D )
term2Pareto <- ddply(grid, .(q), function(data) with(data, price*(mean(pmin(q,D)*(D >= a)))))$V1
term2 <- with(mixture, PLT$w*term2PLT + Pareto$w*term2Pareto)

dataPlot <- rbind(dataPlot, data.frame(q, term1, value = term1 + term2, method = "Mixture of Pareto and 2-PLT"))

# Visualize bound
bitmap("pics/newsVendor_bound_logNorm_70th_percentile.tiff",res = 300, width = 5,height = 5)
library(grid)
library(ggplot2)
plot <- ggplot(dataPlot, aes(x = q, y = value)) + 
  geom_line(aes(linetype = method)) + 
  labs(x = "q", y = "Expected Profit", linetype = "") + 
  geom_vline(xintercept = a, linetype = "dashed") +
  theme(legend.position = c(1/2, 1/8)) 


gtext <-textGrob("a", y = -0.02)
plot <- ggplotGrob(plot + annotation_custom(gtext,xmin = 29, xmax = 29, ymin = -Inf, ymax = Inf))
plot$layout$clip[plot$layout$name=="panel"] <- "off"
grid.draw(plot)
dev.off()
# The command below only works for Mac OS X systems
# It converts to a png format without loss 
# system("sips -s format png pics/newsVendor_bound_logNorm_70th_percentile.tiff --out pics/newsVendor_bound_logNorm_70th_percentile.png") # To convert in a png format without loss 

ggsave(plot,file = "pics/Figure10_newsVendor_bound_logNorm_70th_percentile.svg", width = 5,height = 5,dpi=300)
