library(ggplot2)
library(cowplot)

##correction plots
diagnostic <- read.csv("Documents/FIBR_temporal/cvtk/real_data/correction_diagnostics_uphased.csv", header=T)
head(diagnostic)
table(diagnostic$correction)
table(diagnostic$seqid)
summary(diagnostic)

correction_t <- subset(diagnostic, correction == "True")
correction_f <- subset(diagnostic, correction == "False")



correction_cov_plot <- ggplot(data=correction_t,aes(x=depth, y=offdiag)) +
  geom_point(aes(color=seqid))+
  xlab("Depth")+
  ylab("Covariance")+
  ggtitle("Bias correction")+
  theme_bw()+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position="none")+
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black")

No_correction_cov_plot <- ggplot(data=correction_f,aes(x=depth, y=offdiag)) +
  geom_point(aes(color=seqid))+
  xlab("Depth")+
  ylab("Covariance")+
  theme_bw()+
#  ggtitle("No bias correction")+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position="none")+
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black")

head(correction_f)

No_correction_var_plot <- ggplot(data=correction_f,aes(x=depth, y=diag)) +
  geom_point(aes(color=seqid))+
  xlab("Depth")+
  ylab("Variance")+
  theme_bw()+
#  ggtitle("No bias correction")+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position="none")+
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black")


correction_var_plot <- ggplot(data=correction_t,aes(x=depth, y=diag)) +
  geom_point(aes(color=seqid))+
  xlab("Depth")+
  ylab("Variance")+
  theme_bw()+
  ggtitle("Bias correction")+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position="none")+
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black")


correction_var_plot

pdf("Documents/FIBR_temporal/cvtk/diagnostic_plots.pdf", width=10, height = 6)
plot_grid(No_correction_var_plot, correction_var_plot, No_correction_cov_plot, correction_cov_plot)
dev.off()

final_plot <- plot_grid(No_correction_var_plot, No_correction_cov_plot)
final_plot

png("Documents/FIBR_temporal/cvtk/real_data/diagnostic_plots_final.png", units = "in", width=10, height=5, res=1200)
plot_grid(No_correction_var_plot, No_correction_cov_plot)
dev.off()
