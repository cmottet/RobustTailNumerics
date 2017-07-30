remove(list = ls())
library(svglite)
library(ggplot2)
library(plyr)
library(dplyr)

# Load data
load("data/syntDatalogNormal.RData")
meanlog <- 0
sdlog <- 0.5

###
### Get the data for Figure 12, i.e. the 95% confidence intervals
### for eta, nu, and beta
###
a <- seq(0.8,5, length =100)
CI <- RobustTail::getCIMomentAndDerivatives(sample,
                                     a,
                                     m = 0,
                                     d= 1:2,
                                     nboot = 1E3,
                                     alpha = 0.05,
                                     bootSample = TRUE,
                                     mc.cores = 1) # Can be increased to run computations in parallel
save(CI,file = "data/syntDatalogNormalCI.RData")

###
### Reformat CI as a data frame
###
load("data/syntDatalogNormalCI.RData")
dataPlot <- NULL
for (i in 1:length(CI))
{
  a <- CI[[i]]$a
  truth <- c(DistributionPty::Dlnorm(a,2,meanlog,sdlog),
             DistributionPty::Dlnorm(a,1,meanlog,sdlog),
             1 - DistributionPty::Dlnorm(a,0,meanlog,sdlog))

  newDataPlot <- data.frame(a = a,
                            parameter = rep(c("Density derivative function", "Density function", "Tail distribution function"),3),
                            value = as.numeric(c(CI[[i]]$hyperrectangle[1,], CI[[i]]$hyperrectangle[2,], truth)),
                            group =  rep(c("lB","uB","truth"),each = 3),
                            type = c(rep("Bootstrap 95% CI", 6), rep("True function",3)))
  dataPlot <- rbind(dataPlot, newDataPlot)
}


##
## Plot CI's
##
bitmap("pics/FitKEModelLogNorm.tiff",res = 300, width = 5,height = 5)
plot<-ggplot(dataPlot, aes(x = a, y  = value, group = group)) +
  geom_line(aes(linetype = type)) +
  labs(y = "", linetype = "", x = "") +
  facet_wrap(~parameter, ncol = 2, scales = "free") +
  theme(legend.position = c(7/8, 1/8), legend.justification = c(1, 0))
dev.off()
# The command below only works for Mac OS X systems
# It converts to a png format without loss
# system("sips -s format png pics/FitKEModelLogNorm.tiff --out pics/FitKEModelLogNorm.png")

ggsave(plot,file = "pics/Figure12_FitKEModelLogNorm.svg", width = 5,height = 5,dpi=300)
ggsave(plot,file = "pics/Figure12_FitKEModelLogNorm.pdf", width = 5,height = 5,dpi=300)

