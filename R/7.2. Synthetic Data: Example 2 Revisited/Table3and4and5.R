remove(list=ls())

library(plyr)
library(dplyr)
library(xtable)

###
### Load data
### 

load("data/runSyntheticLogNormalProbCovGPD.RData")
load("data/runSyntheticLogNormalProbCov.RData")

meanlog <- 0   
sdlog <- 0.5

# Function to compute the average upper bound
# and the probability of coverage for data driven approach
getProbCovandAvgUpperBound <- function(table){
  ddply(table, .(a,c), function(data){
    a <- unique(data$a)
    c <- unique(data$c)
    
    truth <- diff(plnorm(c(c,c+1),meanlog,sdlog))
    avgBound <- mean(data$bound, na.rm = TRUE)
    probCov <- mean( diff(plnorm(c(c,c+1),meanlog,sdlog)) <= data$bound, na.rm = TRUE)
    output <- data.frame(d = c+1,truth,avgBound,probCov)
    names(output) <- c("d", "Truth", "Mean upper bound", "Coverage probability")
    return(output)
  })
}

# Create Table 3 in Latex exportable format (a = 3.1)
caption <- "Mean upper bounds and empirical coverage probabilities using worst-case approach with threshold $a = 3.1$."
label <- "Tab:robustnessOptim"
display <- c("d","d","d","E","E","f")
align <- rep("c", length(display))

tabOptimBound31 <- getProbCovandAvgUpperBound(filter(optimBound, a == 3.1))
tabOptimBound31 %>% 
  select(-a) %>%
  xtable(caption = caption, label = label, display = display,align = align) %>%
  print(include.rownames=FALSE, type = "latex",file = "tables/ProbCov1.tex")

# Create Table 4 in Latex exportable format (a = 2.8)
caption <- "Mean upper bounds and empirical coverage probabilities using worst-case approach with threshold $a = 2.8$."
label <- "Tab:robustnessOptim2"

tabOptimBound28 <- getProbCovandAvgUpperBound(filter(optimBound, a == 2.8))
tabOptimBound28 %>% 
  select(-a) %>%
  xtable(caption = caption, label = label,display = display) %>%
  print(include.rownames=FALSE, type = "latex",file = "tables/ProbCov2.tex")

# Create Table 5 in Latex exportable format (GPD approach)
caption <- "Mean upper bounds and empirical coverage probabilities using GPD approach."
label <- "Tab:robustnessGPD"
display <- c("d","d","d","E","E","f")
align <- rep("c", length(display))

tabGPDBound <- getProbCovandAvgUpperBound(transform(GPDBound,a = "u")) %>% select(-a)
tabGPDBound %>%
  xtable(caption = caption, label = label, display = display,align = align) %>%
  print(include.rownames=FALSE, type = "latex",file = "tables/ProbCovGPD.tex")

