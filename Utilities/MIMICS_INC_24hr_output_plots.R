############
# MIMout MIMICS_INC_hourly Plots
###########

library(ggplot2)
library(ggpubr)

# Litter mass
plot_LIT <- ggplot(MIMout, aes(y=LITs, x=DAY, color="Structural")) + geom_line(size=1) +
  geom_line(aes(y=LITm, x=DAY, color="Metabolic"), size=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")

# SOM & MIC pools
plot_SOM_MIC <- ggplot(MIMout, aes(SOMc, x=DAY, color="SOMc")) + geom_line(size=1) +
  geom_line(aes(y=SOMa, x=DAY, color="SOMa"), size=1) +
  geom_line(aes(y=MICr, x=DAY, color="MIC-r"), size=1) +
  geom_line(aes(y=MICK, x=DAY, color="MIC-K"), size=1) +
  theme_bw() +
  ylab("Microbial and soil C") +
  xlab("Incubation Time (days)") +
  labs(color = "C Pool") +
  ylim(0, 3)

# CO2 fraction
plot_CO2 <- ggplot(MIMout, aes(y=rowSums(MIMout[,10:11])/rowSums(MIMout[,3:11]),
                   x=DAY, color="CO2-C")) + geom_line(size=1) +
  theme_bw() +
  ylab("CO2 (fraction of initial") +
  xlab("Incubation Time (days)") +
  labs(color = "C Pool")


# Build a panel plot
ggarrange(plot_LIT, plot_SOM_MIC, plot_CO2,
          nrow=3,
          ncol=1)