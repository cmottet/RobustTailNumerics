remove(list = ls())


load ("data/syntDatalogNormalCI.RData")
load("data/syntDatalogNormal.RData")

# Parameters of the problem
c <- 4
d <- 5
a <- CI[[55]]$a# A index is 55 (compare to the full vector a)  # See Analysis Synthetic Log Normal to get this
H <- function(x) output <- 1/2*( (x+a -c)^2*(x+a>=c) - (x+a -d)^2*(x+a>=d))*(x >=a)

# Case 1 : Data driven optimal bound
etaCI <- CI[[55]]$hyperrectangle$d1
betaCI <- CI[[55]]$hyperrectangle$m0
nuCI <- -CI[[55]]$hyperrectangle$d2[1]
muCI <- etaCI/nuCI
sigmaCI <- 2*betaCI/nuCI
optimBoundCI <- RobustTail::computeBound(H, muCI, sigmaCI, lambda = 0, nuCI)

# Case 2 : True parameters optimal bound
library(DistributionPty) 
meanlog <- 0
sdlog <- 0.5
beta <- 1 - plnorm(a,meanlog,sdlog)
eta <- Dlnorm(a,1,meanlog,sdlog)
nu <- -Dlnorm(a,2, meanlog,sdlog) 
mu <- eta/nu
sigma <- 2*beta/nu
optimBoundVal <- RobustTail::computeBound(H, mu, sigma, lambda = 0,  nu)

# Case 3 : GPD approach 
library(dplyr)
u <- 1.8 
QRM::MEplot(sample) %>% abline(v = u) # Choose the threshold

fitGPD <- QRM::fit.GPD(sample,threshold = u,type = "ml") # Fit the GPD
h <- function(xi, beta) QRM::pGPD(d -u,xi,beta) - QRM::pGPD(c -u,xi,beta)
hGrad <- function(xi,beta) gradientpGPD(d-u,xi,beta) - gradientpGPD(c-u,xi,beta)
GPDbound <- asymptoticCIforGPDfit(fitGPD,h,hGrad,verbose = FALSE) # Compute the CI

# Create Table 1 in Latex exportable format
library(xtable)
Fn <- ecdf(sample)
truth <- diff(plnorm(c(c,d),meanlog,sdlog))

table <- data.frame(c("Truth","ECDF", "GPD","Worst-case with known parameters", "Worst-case appoach"),
                    c(truth, Fn(c) - Fn(d), GPDbound$uB, optimBoundVal$bound, optimBoundCI$bound)) 
names(table)<- c("Method", "Estimated upper bound")

caption <- "Estimated upper bounds of the probability $P(4<X<5)$ for the synthetic data in Example \\ref{example:synthetic}."
label <- "Tab:Ex SyntD"
display <- c("d","s","E")
align  <- rep("c",length(display))
xtable(table, caption = caption , label = label, display = display, align = align) %>% 
  print(include.rownames=FALSE, type = "latex",file = "tables/logNormal45.tex")
