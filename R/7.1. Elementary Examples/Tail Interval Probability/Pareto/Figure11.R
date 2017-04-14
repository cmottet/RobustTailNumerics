remove(list= ls())

# Create the data for the plot
library(DistributionPty)
library(dplyr)
library(ggplot2)
library(svglite)

shape = scale = 1
load("data/runSyntheticPareto.RData")
widthRect <- diff(ppareto(optimBound$a, shape,scale)[1:2])
dataPlot <- transform(optimBound,
                      aLperc = ppareto(a, shape,scale),
                      aUperc = ppareto(a, shape,scale) + widthRect,
                      pc = ppareto(c, shape,scale),
                      pd = ppareto(d, shape,scale),
                      truth =  ppareto(d, shape,scale) - ppareto(c, shape,scale)) %>%
             transform(ratio = bound/truth)

# Create the Plot
bitmap("pics/runSyntheticPareto.tiff",res = 300, width = 5,height = 5)
plot <- ggplot(dataPlot, aes(xmin = aLperc, xmax = aUperc,  ymin = pc, ymax = pd, fill = ratio)) +
  geom_rect(colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  scale_x_continuous(breaks = seq(0.7,0.86, length = 9),expand=c(0,0)) +
  scale_y_continuous(breaks = seq(0.86,0.98, length = 7),expand=c(0,0)) +
  labs(y = "Interval (c,d) in percentile",
       x = "Threshold a in percentile",
       fill = "Ratio between the optimal upper bound and the true value")  +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 20, title.position = "right"))+
  theme(legend.title = element_text(angle = 90))
ggsave(plot,file= "pics/Figure11_runSyntheticPareto.svg",dpi = 300, width = 5,height = 5)

dev.off()
# The command below only works for Mac OS X systems
# It converts to a png format without loss
# system("sips -s format png pics/runSyntheticPareto.tiff --out pics/runSyntheticPareto.png")
