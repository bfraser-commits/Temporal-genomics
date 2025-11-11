library(ggplot2)
library(DataCombine)
library(cowplot)

test <- read.csv("Documents/Mijke/FIBR_Ms/revisions/genomicscensusdata.csv", header=T)


total_f <- data.frame(stream = test$stream, sampling = test$sampling, N.hat = test$F.N.hat, N.Se = test$F.N.se, Sex = "F")
total_m <- data.frame(stream = test$stream, sampling = test$sampling, N.hat = test$M.N.hat, N.Se = test$M.N.se, Sex = "M")


test$sum <- test$M.N.hat + test$F.N.hat
testnew <- subset(test, sampling < 120)

sum_text <- data.frame(stream= c("CA", "TA", "LL", "UL"), x=c(121,121,121,121), y=c(367,893,1200,1613))



summary(testnew$sampling)

#scale_x_continuous(limits = c(0,120),breaks = c(1,23,47,71,95,119), labels=c(2008,2010,2012,2014,2016,2018))+


###plot with year on top and month on bottom
subset(testnew, sampling == 29)

regulation_dots <- data.frame(stream= c("CA", "TA", "LL", "UL"), x=c(21,15,32,29), y=c(213.3311,274.9192,1083.5503,2612.1831))

two_axis_plot <- ggplot(data = testnew, aes(x = sampling, y = sum)) +
  geom_line(aes(color = stream), size = 1.5, alpha = 0.8) +
  geom_text(data = sum_text,
            mapping = aes(x = x, y = y, label = stream),
            size = 6,
            colour = c(CA = "#FFB90F", TA = "#FF7F00", LL = "#5D478B", UL = "#CD0000")) +
  # Annotations
  annotate('segment', x = 35, xend = 35, y = 0, yend = 3000, size = 0.5, linetype = "dashed")+
  annotate("segment", x = 59, xend = 59, y = 0, yend = 3200, size = 0.5, linetype = "dashed") +
  annotate("segment", x = 119, xend = 119, y = 0, yend = 3200, size = 0.5, linetype = "dashed") +
  annotate('text', x= 35, y = 3200, label = 'male size change', size = 4.5, hjust = 0.5)+
  annotate("text", x = 59, y = 3200, label = "genome data", size = 4.5, hjust = 0.5) +
  annotate("text", x = 59, y = 3000, label = "male colour change", size = 4.5, hjust = 0.5) +
  annotate("text", x = 119, y = 3200, label = "genome data", size = 4.5, hjust = 0.5) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.position = "none",
    plot.margin = margin(c(0.5, 0.7, 0.1, 0.1), unit = "cm")
  ) +
  
  scale_colour_manual(values = c(CA = "#FFB90F", TY = "#FF7F00", LL = "#5D478B", UL = "#CD0000")) +
  
  xlab(expression("Sampling Month")) +
  ylab(expression("N (total)")) +
  
  # Add top x-axis: convert sampling months to years
  scale_x_continuous(
    limits = c(0, 122),
    breaks = c(1, 20, 40, 60, 80, 100, 120),
    labels = c(1, 20, 40, 60, 80, 100, 120),
    sec.axis = dup_axis(  # top axis
      trans = ~ .,
      breaks = c(1, 23, 47, 71, 95, 119),
      labels = c(2008, 2010, 2012, 2014, 2016, 2018),
      name = "Year"
    )
  ) +
  
  scale_y_continuous(limits = c(0, 3200))

pdf("Documents/FIBR_temporal/census_plot_annot.pdf", width = 10, height=5)
two_axis_plot
dev.off()
