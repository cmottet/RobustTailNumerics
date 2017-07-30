remove(list = ls())
library(svglite)
library(plyr)
library(dplyr)
library(ggplot2)

# Load Danish insurance fire data
load("data/datafireinsurance.RData")
sample <- as.numeric(sample$Loss.in.DKM)

amax <- sort(sample)[length(sample) -15]

###
### Compute the 95% confidence intervals
### for the parameters eta, beta, nu over different
### values of a
###
a <- seq(1, amax, length = 50)
CI <- RobustTail::getCIMomentAndDerivatives(sample,
                                     a,
                                     m = 0,
                                     d= 1:2,
                                     nboot = 1E3,
                                     alpha = 0.05,
                                     bootSample = TRUE,
                                     mc.cores = 1) # Increase the number of cores to use parallel computing
#save(CI,file = "data/datafireinsuranceCI.RData")
load("data/datafireinsuranceCI.RData")
###
### Transform the CI's in a data frame format
###

dataPlot <- NULL
for (i in 1:length(CI))
{
  bootSample <- CI[[i]]$bootSample

  newDataPlot <- data.frame(a =  CI[[i]]$a,
                            parameter = rep(c("Density derivative function", "Density function", "Tail distribution function"),3),
                            value =  as.numeric(c(CI[[i]]$hyperrectangle[1,], CI[[i]]$hyperrectangle[2,], as.numeric(apply(bootSample,2,mean)))),
                            group = rep(c("lB", "uB", "Fhat"),each  = 3),
                            type = c(rep("Bootstrap 95% CI", 6),rep("Bootstrap-estimated function", 3) ))
  dataPlot <- rbind(dataPlot, newDataPlot)
}

###
### Plot
###
bitmap("pics/FitKEFire.tiff",res = 300, width = 5,height = 5)
plot<-ggplot(dataPlot, aes(x = a, y  = value, group = group)) +
  geom_line(aes(linetype = type)) +
  labs(y = "", linetype = "", x = "") +
  facet_wrap(~parameter, ncol = 2, scales = "free") +
  theme(legend.position = c(7/8, 1/8), legend.justification = c(0.75, 0))
dev.off()
# The command below only works for Mac OS X systems
# It converts to a png format without loss
# system("sips -s format png pics/FitKEFire.tiff --out pics/FitKEFire.png")

ggsave(plot,file = "pics/Figure13_FitKEFire.svg", width = 5,height = 5,dpi=300)
ggsave(plot,file = "pics/Figure13_FitKEFire.pdf", width = 5,height = 5,dpi=300)

