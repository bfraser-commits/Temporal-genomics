library(ggplot2)
library(cowplot)
#SLIMS
real <- read.table("Documents/FIBR_temporal/cvtk/real_data/gw_covs_mean_unphased.tsv", header=T,  sep ="\t")

CA_time <-real[2, 1]
TA_time <- real[4, 3]
LL_time <- real[6, 5]
UL_time <- real[8, 7]

CA_time
CA_time_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_time.txt")
CA_time_percentile <- mean(CA_time_slims$V2 <= CA_time)
CA_time_percentile


TA_time_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_time.txt")
hist(TA_time_slims$V2)
TA_time_percentile <- mean(TA_time_slims$V2 <= TA_time)
TA_time_percentile

LL_time_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/LL_time.txt")
hist(LL_time_slims$V2)
LL_time_percentile <- mean(LL_time_slims$V2 <= LL_time)
LL_time_percentile

UL_time_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/UL_time.txt")
hist(UL_time_slims$V2)
UL_time_percentile <- mean(UL_time_slims$V2 <= UL_time)
UL_time_percentile


##plot
all_time <- data.frame("covariance" = c(CA_time, TA_time, LL_time, UL_time), "percentile" = c(CA_time_percentile, TA_time_percentile, LL_time_percentile, UL_time_percentile), "pop" = c("CA", "TA", "LL", "UL"), "slims" = "NO")
all_time
CA_time_slims_df <- data.frame("covariance" = CA_time_slims$V2, "pop" = "CA", "slims" = "SLIM")
TA_time_slims_df <- data.frame("covariance" = TA_time_slims$V2, "pop" = "TA", "slims" = "SLIM")
LL_time_slims_df <- data.frame("covariance" = LL_time_slims$V2, "pop" = "LL", "slims" = "SLIM")
UL_time_slims_df <- data.frame("covariance" = UL_time_slims$V2, "pop" = "UL", "slims" = "SLIM")
all_time_slims <- data.frame(rbind(CA_time_slims_df, TA_time_slims_df, LL_time_slims_df, UL_time_slims_df))
summary(all_time_slims)
all_time_slims$pop <- factor(all_time_slims$pop, levels = c("CA", "TA", "LL", "UL"))
all_time$pop <- factor(all_time$pop, levels = c("CA", "TA", "LL", "UL"))
summary(all_time_slims)

all_time_plot <- ggplot() +
  geom_violin(
    data = all_time_slims,
    aes(x = pop, y = covariance),
    trim = FALSE, width = 0.4
  ) +
  geom_boxplot(
    data = all_time_slims,
    aes(x = pop, y = covariance),
    width = 0.1, alpha = 0.8, outlier.shape = NA,
    color = "black", fill = NA
  ) +
  geom_point(
    data = all_time,
    aes(
      x = as.numeric(factor(pop)) + 0.3,
      y = covariance,
      fill = percentile
    ),
    shape = 21,          
    colour = "black",    
    stroke = 0.8,         
    size = 3,
    alpha = 0.95
  ) +
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.2))) +
  labs(title = "Across Time", x = "Population", y = "Covariance") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16)
  )


pdf("Documents/FIBR_temporal/cvtk/across_time_cov_violin_colour.pdf", width=6, height=4)
all_time_plot
dev.off()

##time 1
CA_TA_13 <- real[3, 1]
CA_TA_13
CA_LL_13 <- real[5, 1]
CA_LL_13
CA_UL_13 <- real[7, 1]
CA_UL_13
TA_LL_13 <- real[5, 3]
TA_LL_13
TA_UL_13 <- real[7, 3]
TA_UL_13
LL_UL_13 <- real[7, 5]
LL_UL_13

CA_TA_13_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_TA_1.txt")
hist(CA_TA_13_slims$V2)
CA_TA_13_percentile <- mean(CA_TA_13_slims$V2 <= CA_TA_13)
CA_TA_13_percentile

CA_LL_13_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_LL_1.txt")
hist(CA_LL_13_slims$V2)
CA_LL_13_percentile <- mean(CA_LL_13_slims$V2 <= CA_LL_13)
CA_LL_13_percentile

CA_UL_13_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_UL_1.txt")
hist(CA_UL_13_slims$V2)
CA_UL_13_percentile <- mean(CA_UL_13_slims$V2 <= CA_UL_13)
CA_UL_13_percentile

TA_LL_13_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_LL_1.txt")
hist(TA_LL_13_slims$V2)
TA_LL_13_percentile <- mean(TA_LL_13_slims$V2 <= TA_LL_13)
TA_LL_13_percentile

TA_UL_13_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_UL_1.txt")
hist(TA_UL_13_slims$V2)
TA_UL_13_percentile <- mean(TA_UL_13_slims$V2 <= TA_UL_13)
TA_UL_13_percentile

LL_UL_13_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/LL_UL_1.txt")
hist(LL_UL_13_slims$V2)
LL_UL_13_percentile <- mean(LL_UL_13_slims$V2 <= LL_UL_13)
LL_UL_13_percentile



##plot
all_13 <- data.frame("covariance" = c(CA_TA_13, CA_LL_13, CA_UL_13, TA_LL_13, TA_UL_13, LL_UL_13), "percentile" = c(CA_TA_13_percentile, CA_LL_13_percentile, CA_UL_13_percentile, TA_LL_13_percentile, TA_UL_13_percentile, LL_UL_13_percentile), "comp" = c("CA_TA", "CA_LL", "CA_UL", "TA_LL", "TA_UL", "LL_UL"), "slims" = "NO")
all_13
CA_TA_13_slims_df <- data.frame("covariance" = CA_TA_13_slims$V2, "comp" = "CA_TA", "slims" = "SLIM")
CA_LL_13_slims_df <- data.frame("covariance" = CA_LL_13_slims$V2, "comp" = "CA_LL", "slims" = "SLIM")
CA_UL_13_slims_df <- data.frame("covariance" = CA_UL_13_slims$V2, "comp" = "CA_UL", "slims" = "SLIM")
TA_LL_13_slims_df <- data.frame("covariance" = TA_LL_13_slims$V2, "comp" = "TA_LL", "slims" = "SLIM")
TA_UL_13_slims_df <- data.frame("covariance" = TA_UL_13_slims$V2, "comp" = "TA_UL", "slims" = "SLIM")
LL_UL_13_slims_df <- data.frame("covariance" = LL_UL_13_slims$V2, "comp" = "LL_UL", "slims" = "SLIM")

all_13_slims <- data.frame(rbind(CA_TA_13_slims_df, CA_LL_13_slims_df, CA_UL_13_slims_df, TA_LL_13_slims_df, TA_UL_13_slims_df, LL_UL_13_slims_df))
summary(all_13_slims)
all_13_slims$comp <- factor(all_13_slims$comp, levels = c("CA_TA", "CA_LL", "CA_UL", "TA_LL", "TA_UL", "LL_UL"))
all_13$comp <- factor(all_13$comp, levels = c("CA_TA", "CA_LL", "CA_UL", "TA_LL", "TA_UL", "LL_UL"))



all_13_plot <- ggplot() +
  geom_violin(data = all_13_slims, 
              aes(x = comp, y = covariance),
              trim = FALSE, width = 0.4) +
  # Boxplot for quartiles - narrower
  geom_boxplot(data = all_13_slims,
               aes(x = comp, y = covariance),
               width = 0.1, alpha = 0.8, outlier.shape = NA,
               color = "black", fill = NA) +
  # Observed points - further to the right
  geom_point(
    data = all_13,
    aes(
      x = as.numeric(factor(comp)) + 0.3,
      y = covariance,
      fill = percentile
    ),
    shape = 21,          
    colour = "black",    
    stroke = 0.8,         
    size = 3,
    alpha = 0.95
  ) +
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  scale_x_discrete(expand = expansion(add = c(0.6, 1.2))) +
  labs(title = "1st Time Interval",
       x = "Comparison",
       y = "Covariance") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_blank(),
    plot.title = element_text(size = 16)
  )

all_13_plot
summary(all_13$percentile)

pdf("Documents/FIBR_temporal/cvtk/time_interval_1_violin_colour.pdf", width=8, height=4)
all_13_plot
dev.off()

#timepoint2
CA_TA_18 <- real[4, 2]
CA_LL_18 <- real[6, 2]
CA_UL_18 <- real[8, 2]
TA_LL_18 <- real[6, 4]
TA_UL_18 <- real[8, 4]
LL_UL_18 <- real[8, 6]



###
CA_TA_18_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_TA_2.txt")
hist(CA_TA_18_slims$V2)
CA_TA_18_percentile <- mean(CA_TA_18_slims$V2 <= CA_TA_18)
CA_TA_18_percentile

CA_LL_18_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_LL_2.txt")
hist(CA_LL_18_slims$V2)
CA_LL_18_percentile <- mean(CA_LL_18_slims$V2 <= CA_LL_18)
CA_LL_18_percentile

CA_UL_18_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_UL_2.txt")
hist(CA_UL_18_slims$V2)
CA_UL_18_percentile <- mean(CA_UL_18_slims$V2 <= CA_UL_18)
CA_UL_18_percentile

TA_LL_18_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_LL_2.txt")
hist(TA_LL_18_slims$V2)
TA_LL_18_percentile <- mean(TA_LL_18_slims$V2 <= TA_LL_18)
TA_LL_18_percentile

TA_UL_18_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_UL_2.txt")
hist(TA_UL_18_slims$V2)
TA_UL_18_percentile <- mean(TA_UL_18_slims$V2 <= TA_UL_18)
TA_UL_18_percentile

LL_UL_18_slims <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/LL_UL_2.txt")
hist(LL_UL_18_slims$V2)
LL_UL_18_percentile <- mean(LL_UL_18_slims$V2 <= LL_UL_18)
LL_UL_18_percentile

#plot
all_18 <- data.frame("covariance" = c(CA_TA_18, CA_LL_18, CA_UL_18, TA_LL_18, TA_UL_18, LL_UL_18), "percentile" = c(CA_TA_18_percentile, CA_LL_18_percentile, CA_UL_18_percentile, TA_LL_18_percentile, TA_UL_18_percentile, LL_UL_18_percentile), "comp" = c("CA_TA", "CA_LL", "CA_UL", "TA_LL", "TA_UL", "LL_UL"), "slims" = "NO")
all_18
CA_TA_18_slims_df <- data.frame("covariance" = CA_TA_18_slims$V2, "comp" = "CA_TA", "slims" = "SLIM")
CA_LL_18_slims_df <- data.frame("covariance" = CA_LL_18_slims$V2, "comp" = "CA_LL", "slims" = "SLIM")
CA_UL_18_slims_df <- data.frame("covariance" = CA_UL_18_slims$V2, "comp" = "CA_UL", "slims" = "SLIM")
TA_LL_18_slims_df <- data.frame("covariance" = TA_LL_18_slims$V2, "comp" = "TA_LL", "slims" = "SLIM")
TA_UL_18_slims_df <- data.frame("covariance" = TA_UL_18_slims$V2, "comp" = "TA_UL", "slims" = "SLIM")
LL_UL_18_slims_df <- data.frame("covariance" = LL_UL_18_slims$V2, "comp" = "LL_UL", "slims" = "SLIM")

all_18_slims <- data.frame(rbind(CA_TA_18_slims_df, CA_LL_18_slims_df, CA_UL_18_slims_df, TA_LL_18_slims_df, TA_UL_18_slims_df, LL_UL_18_slims_df))
summary(all_18_slims)
all_18_slims$comp <- factor(all_18_slims$comp, levels = c("CA_TA", "CA_LL", "CA_UL", "TA_LL", "TA_UL", "LL_UL"))
all_18$comp <- factor(all_18$comp, levels = c("CA_TA", "CA_LL", "CA_UL", "TA_LL", "TA_UL", "LL_UL"))



all_18_plot <- ggplot() +
  geom_violin(data = all_18_slims, 
              aes(x = comp, y = covariance),
              trim = FALSE, width = 0.4) +
  geom_boxplot(data = all_18_slims,
               aes(x = comp, y = covariance),
               width = 0.1, alpha = 0.8, outlier.shape = NA,
  ) +
  geom_point(
    data = all_18,
    aes(
      x = as.numeric(factor(comp)) + 0.3,
      y = covariance,
      fill = percentile
    ),
    shape = 21,          
    colour = "black",    
    stroke = 0.8,         
    size = 3,
    alpha = 0.95
  ) +
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )  +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.2))) +
  labs(title = "2nd Time Interval",
       x = "Comparison",
       y = "Covariance") +
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_blank(),
    plot.title = element_text(size = 16)
  )
all_18_plot

pdf("Documents/FIBR_temporal/cvtk/time_interval_2_violin_colour.pdf", width=8, height=4)
all_18_plot
dev.off()
