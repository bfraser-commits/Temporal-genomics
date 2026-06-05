library(ggplot2)
library(cowplot)

#CA_time
CA_time <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_time_chr.txt")
CA_time$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_time$V1)
CA_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_time.txt")
CA_slim$chr <- ("SLIM")
CA_time
CA_time$chr <- factor(CA_time$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

head(CA_slim)
head(CA_time)
CA_time$percentile <- sapply(CA_time$V2, function(v) mean(CA_slim$V2 <= v))

# plot
CA_time_plot <- ggplot() +
  geom_violin(data = CA_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = CA_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_time,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
    scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "CA") +
  ylim(-0.035,0.06) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_time_plot
subset(CA_time, percentile < 0.05)

#TA_time
TA_time <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/TA_time_chr.txt")
TA_time$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", TA_time$V1)
TA_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_time.txt")
TA_slim$chr <- ("SLIM")
TA_time$chr <- factor(TA_time$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

TA_time$percentile <- sapply(TA_time$V2, function(v) mean(TA_slim$V2 <= v))
TA_time
subset(TA_time, percentile < 0.05)

# plot
TA_time_plot <- ggplot() +
  geom_violin(data = TA_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = TA_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = TA_time,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "TA") +
  ylim(-0.035,0.06) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

TA_time_plot



#LL
LL_time <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/LL_time_chr.txt")
LL_time$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", LL_time$V1)
head(LL_time)
LL_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/LL_time.txt")
LL_slim$chr <- ("SLIM")
LL_time$chr <- factor(LL_time$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

LL_time$percentile <- sapply(LL_time$V2, function(v) mean(LL_slim$V2 <= v))
subset(LL_time, percentile < 0.05)

# plot
LL_time_plot <- ggplot() +
  geom_violin(data = LL_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = CA_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = LL_time,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "LL") +
  ylim(-0.035,0.06) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

LL_time_plot

nrow(subset(LL_time, LL_time$percentile > 0.95))

#UL
UL_time <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/UL_time_chr.txt")
UL_time$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", UL_time$V1)
head(UL_time)
UL_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/UL_time.txt")
UL_slim$chr <- ("SLIM")
UL_time$chr <- factor(UL_time$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

UL_time$percentile <- sapply(UL_time$V2, function(v) mean(UL_slim$V2 <= v))
subset(UL_time, percentile < 0.05)

# plot
UL_time_plot <- ggplot() +
  geom_violin(data = UL_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = UL_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = UL_time,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "UL") +
  ylim(-0.035,0.06) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

UL_time_plot
#subset(UL_time, UL_time$percentile > 0.95)
#subset(TA_time, TA_time$percentile > 0.95)


png("Documents/FIBR_temporal/cvtk/across_time_cov_violin_chrs_colour.png", width=10, height=6, units = "in", res = 300)
plot_grid(CA_time_plot, TA_time_plot, LL_time_plot, UL_time_plot)
dev.off()

###time point 1##
#CA_TA_1
CA_TA_1 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_TA_1_chr.txt")
CA_TA_1$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_TA_1$V1)
CA_TA_1_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_TA_1.txt")
CA_TA_1_slim$chr <- ("SLIM")
CA_TA_1$chr <- factor(CA_TA_1$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))


CA_TA_1$percentile <- sapply(CA_TA_1$V2, function(v) mean(CA_TA_1_slim$V2 <= v))

# plot
CA_TA_1_plot <- ggplot() +
  geom_violin(data = CA_TA_1_slim,
              aes(x = chr, y = V2))+
  geom_boxplot(data = CA_TA_1_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_TA_1,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "CA_TA_1") +
  ylim(-0.009,0.09) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_TA_1_plot
subset(CA_TA_1, CA_TA_1$percentile > 0.95)
nrow(subset(CA_TA_1, CA_TA_1$percentile > 0.95))

#CA_TA_1
CA_LL_1 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_LL_1_chr.txt")
CA_LL_1$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_LL_1$V1)
CA_LL_1_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_LL_1.txt")
CA_LL_1_slim$chr <- ("SLIM")
CA_LL_1$chr <- factor(CA_LL_1$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

CA_LL_1$percentile <- sapply(CA_LL_1$V2, function(v) mean(CA_LL_1_slim$V2 <= v))

# plot
CA_LL_1_plot <- ggplot() +
  geom_violin(data = CA_LL_1_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = CA_LL_1_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_LL_1,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "CA_LL_1") +
  ylim(-0.009,0.09)+
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_LL_1_plot
subset(CA_LL_1, CA_LL_1$percentile > 0.90)
nrow(subset(CA_LL_1, CA_LL_1$percentile > 0.95))

CA_LL_1
#CA_LL_1$V2[CA_LL_1$chr == "chr15"]
#mean(CA_LL_1_slim$V2 <= CA_LL_1$V2[CA_LL_1$chr == "chr15"])


CA_UL_1 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_UL_1_chr.txt")
CA_UL_1$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_UL_1$V1)
CA_UL_1_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_UL_1.txt")
CA_UL_1_slim$chr <- ("SLIM")
CA_UL_1$chr <- factor(CA_UL_1$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

CA_UL_1$percentile <- sapply(CA_UL_1$V2, function(v) mean(CA_UL_1_slim$V2 <= v))
CA_UL_1

# plot
CA_UL_1_plot <- ggplot() +
  geom_violin(data = CA_UL_1_slim,
              aes(x = chr, y = V2))+
  geom_boxplot(data = CA_UL_1_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_UL_1,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "CA_UL_1") +
    ylim(-0.009,0.09)+
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_UL_1_plot
subset(CA_UL_1, CA_UL_1$percentile > 0.95)

CA_UL_1
#mean(CA_UL_1_slim$V2 <= CA_UL_1$V2[CA_UL_1$chr == "chr15"])

TA_LL_1 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/TA_LL_1_chr.txt")
TA_LL_1$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", TA_LL_1$V1)
TA_LL_1_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_LL_1.txt")
TA_LL_1_slim$chr <- ("SLIM")
TA_LL_1$chr <- factor(TA_LL_1$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

TA_LL_1$percentile <- sapply(TA_LL_1$V2, function(v) mean(TA_LL_1_slim$V2 <= v))

# plot
TA_LL_1_plot <- ggplot() +
  geom_violin(data = TA_LL_1_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = TA_LL_1_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = TA_LL_1,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "TA_LL_1") +
  ylim(-0.009,0.09)+
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

TA_LL_1_plot
#mean(TA_LL_1_slim$V2 <= TA_LL_1$V2[TA_LL_1$chr == "chr15"])
subset(TA_LL_1, TA_LL_1$percentile > 0.95)

TA_UL_1 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/TA_UL_1_chr.txt")
TA_UL_1$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", TA_UL_1$V1)
TA_UL_1_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_UL_1.txt")
TA_UL_1_slim$chr <- ("SLIM")
TA_UL_1$chr <- factor(TA_UL_1$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

TA_UL_1$percentile <- sapply(TA_UL_1$V2, function(v) mean(TA_UL_1_slim$V2 <= v))

# plot
TA_UL_1_plot <- ggplot() +
  geom_violin(data = TA_UL_1_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = TA_UL_1_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = TA_UL_1,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  )+
  labs(x = "Chr", y = "Covariance", title = "TA_UL_1") +
  ylim(-0.009,0.09) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

TA_UL_1_plot
#mean(TA_UL_1_slim$V2 <= TA_UL_1$V2[TA_UL_1$chr == "chr15"])
subset(TA_UL_1, TA_UL_1$percentile > 0.95)

LL_UL_1 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/LL_UL_1_chr.txt")
LL_UL_1$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", LL_UL_1$V1)
LL_UL_1_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/LL_UL_1.txt")
LL_UL_1_slim$chr <- ("SLIM")
LL_UL_1$chr <- factor(LL_UL_1$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

LL_UL_1$percentile <- sapply(LL_UL_1$V2, function(v) mean(LL_UL_1_slim$V2 <= v))

# plot
LL_UL_1_plot <- ggplot() +
  geom_violin(data = LL_UL_1_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = LL_UL_1_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = LL_UL_1,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "LL_UL_1") +
  ylim(-0.009,0.09) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

LL_UL_1_plot
#mean(LL_UL_1_slim$V2 <= LL_UL_1$V2[LL_UL_1$chr == "chr15"])
subset(LL_UL_1, LL_UL_1$percentile > 0.95)


png("Documents/FIBR_temporal/cvtk/time_1_cov_violin_chrs_colours.png", width=12, height=6, units = "in", res = 300)
plot_grid(CA_TA_1_plot, CA_LL_1_plot, CA_UL_1_plot, TA_LL_1_plot, TA_UL_1_plot, LL_UL_1_plot)
dev.off()


dfs <- list(ds1 = CA_TA_1, ds2 = CA_LL_1, ds3 = CA_UL_1, ds4 = TA_LL_1, ds5 = TA_UL_1, ds6 = LL_UL_1)
thr <- 0.95
head(CA_TA_1)
# for each dataset: unique chrs that exceed the threshold
chr_sets <- lapply(dfs, function(d) unique(d$chr[d$percentile > thr]))

# how many datasets each chr appears in (as > thr)
shared_counts <- sort(table(unlist(chr_sets)), decreasing = TRUE)
shared_counts

##timepoint2
#CA_TA_2
CA_TA_2 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_TA_2_chr.txt")
CA_TA_2$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_TA_2$V1)
CA_TA_2_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_TA_2.txt")
CA_TA_2_slim$chr <- ("SLIM")
CA_TA_2$chr <- factor(CA_TA_2$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

#min(CA_UL_2_slim$V2)
#min(CA_TA_2_slim$V2)

CA_TA_2$percentile <- sapply(CA_TA_2$V2, function(v) mean(CA_TA_2_slim$V2 <= v))

# plot
CA_TA_2_plot <- ggplot() +
  geom_violin(data = CA_TA_2_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = CA_TA_2_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_TA_2,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "CA_TA_2") +
  ylim(-0.021, 0.027) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_TA_2_plot
subset(CA_TA_2, CA_TA_2$percentile > 0.95)

CA_TA_2

CA_LL_2 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_LL_2_chr.txt")
CA_LL_2$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_LL_2$V1)
CA_LL_2_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_LL_2.txt")
CA_LL_2_slim$chr <- ("SLIM")
CA_LL_2$chr <- factor(CA_LL_2$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

CA_LL_2$percentile <- sapply(CA_LL_2$V2, function(v) mean(CA_LL_2_slim$V2 <= v))
CA_LL_2

# plot
CA_LL_2_plot <- ggplot() +
  geom_violin(data = CA_LL_2_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = CA_LL_2_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_LL_2,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "CA_LL_2") +
  ylim(-0.021, 0.027) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_LL_2_plot
summary(CA_LL_2)
subset(CA_LL_2, CA_LL_2$percentile > 0.95)

CA_LL_2

CA_UL_2 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/CA_UL_2_chr.txt")
CA_UL_2$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", CA_UL_2$V1)
CA_UL_2_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/CA_UL_2.txt")
CA_UL_2_slim$chr <- ("SLIM")
CA_UL_2$chr <- factor(CA_UL_2$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

CA_UL_2$percentile <- sapply(CA_UL_2$V2, function(v) mean(CA_UL_2_slim$V2 <= v))

# plot
CA_UL_2_plot <- ggplot() +
  geom_violin(data = CA_UL_2_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = CA_UL_2_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = CA_UL_2,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "CA_UL_2") +
  ylim(-0.021, 0.027) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

CA_UL_2_plot
subset(CA_UL_2, CA_UL_2$percentile > 0.95)

TA_LL_2 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/TA_LL_2_chr.txt")
TA_LL_2$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", TA_LL_2$V1)
TA_LL_2_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_LL_2.txt")
TA_LL_2_slim$chr <- ("SLIM")
TA_LL_2$chr <- factor(TA_LL_2$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

TA_LL_2$percentile <- sapply(TA_LL_2$V2, function(v) mean(TA_LL_2_slim$V2 <= v))

# plot
TA_LL_2_plot <- ggplot() +
  geom_violin(data = TA_LL_2_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = TA_LL_2_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = TA_LL_2,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "TA_LL_2") +
  ylim(-0.021, 0.027) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

TA_LL_2_plot
nrow(subset(TA_LL_2, TA_LL_2$percentile > 0.95))

TA_UL_2 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/TA_UL_2_chr.txt")
TA_UL_2$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", TA_UL_2$V1)
TA_UL_2_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/TA_UL_2.txt")
TA_UL_2_slim$chr <- ("SLIM")
TA_UL_2$chr <- factor(TA_UL_2$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))
TA_UL_2$percentile <- sapply(TA_UL_2$V2, function(v) mean(TA_UL_2_slim$V2 <= v))

# plot
TA_UL_2_plot <- ggplot() +
  geom_violin(data = TA_UL_2_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = TA_UL_2_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = TA_UL_2,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "TA_UL_2") +
  ylim(-0.021, 0.027) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

TA_UL_2_plot
subset(TA_UL_2, TA_UL_2$percentile > 0.95)

LL_UL_2 <- read.table("Documents/FIBR_temporal/cvtk/real_data/chrs/LL_UL_2_chr.txt")
LL_UL_2$chr <- sub(".*/(chr[0-9]+)_.*", "\\1", LL_UL_2$V1)
LL_UL_2_slim <- read.table("Documents/FIBR_temporal/cvtk/slims_out_summary_multpops/LL_UL_2.txt")
LL_UL_2_slim$chr <- ("SLIM")
LL_UL_2$chr <- factor(LL_UL_2$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                              "chr6", "chr7", "chr8", "chr9", "chr10", 
                                              "chr11", "chr12", "chr13", "chr14", "chr15", 
                                              "chr16", "chr17", "chr18", "chr19", "chr20",
                                              "chr21", "chr22", "chr23"))

LL_UL_2$percentile <- sapply(LL_UL_2$V2, function(v) mean(LL_UL_2_slim$V2 <= v))

# plot
LL_UL_2_plot <- ggplot() +
  geom_violin(data = LL_UL_2_slim,
              aes(x = chr, y = V2)) +
  geom_boxplot(data = LL_UL_2_slim,
               aes(x = chr, y = V2),
               width = 0.2, outlier.shape = NA, fill = NA) +
  geom_point(data = LL_UL_2,
             aes(x = chr, y = V2, fill= percentile),
             size = 2, alpha=0.95, shape=21)+
  scale_fill_gradientn(
    colours = c("#2b2b2b", "#9e9e9e", "#fdae61", "#d7191c"),
    values  = c(0, 0.6, 0.9, 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.1),
    labels  = function(x) sprintf("%d", round(100 * x)),
    name    = "Significance (0–100)"
  ) +
  labs(x = "Chr", y = "Covariance", title = "LL_UL_2") +
  ylim(-0.021, 0.027) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

LL_UL_2_plot
subset(LL_UL_2, LL_UL_2$percentile > 0.95)

png("Documents/FIBR_temporal/cvtk/time_2_cov_violin_chrs_colours.png", width=12, height=6, units = "in", res = 300)
plot_grid(CA_TA_2_plot, CA_LL_2_plot, CA_UL_2_plot, TA_LL_2_plot, TA_UL_2_plot, LL_UL_2_plot)
dev.off()

dfs <- list(ds1 = CA_TA_2, ds2 = CA_LL_2, ds3 = CA_UL_2, ds4 = TA_LL_2, ds5 = TA_UL_2, ds6 = LL_UL_2)
thr <- 0.95
# for each dataset: unique chrs that exceed the threshold
chr_sets <- lapply(dfs, function(d) unique(d$chr[d$percentile > thr]))
# how many datasets each chr appears in (as > thr)
shared_counts <- sort(table(unlist(chr_sets)), decreasing = TRUE)
shared_counts



#sup data correlate by size
chr_size <- read.table("Documents/FIBR_temporal/cvtk/slim5_slims/STAR.chromosomes.release.fasta.fai")
head(chr_size)
colnames(chr_size) <- c("chr", "size", "cum_size", "60", "61")

test <- merge(merge(merge(chr_size, CA_time, by = "chr"), CA_TA_1, by = "chr"), CA_TA_2, by = "chr")


names(test)[names(test) == "V2.x"] <- "V2_CAtime"
names(test)[names(test) == "V2.y"] <- "V2_CATA1"
names(test)[names(test) == "V2"] <- "V2_CATA2"
head(test)

CA_time_size <- ggplot(test, aes(x = size, y = V2_CAtime)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "chromosome size", y = "Covariance", title="CA") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 10)
  )


CA_TA_1_size <- ggplot(test, aes(x = size, y = V2_CATA1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "chromosome size", y = "Covariance", title="CA_TA_1") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 10)
  )


CA_TA_2_size <- ggplot(test, aes(x = size, y = V2_CATA2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "chromosome size", y = "Covariance", title="CA_TA_2") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 10))

CA_TA_2_size

plot_grid(CA_time_size, CA_TA_1_size, CA_TA_2_size, nrow=1)


L_test <- merge(merge(merge(chr_size, LL_time, by = "chr"), LL_UL_1, by = "chr"), LL_UL_2, by = "chr")
names(L_test)[names(L_test) == "V2.x"] <- "V2_LLtime"
names(L_test)[names(L_test) == "V2.y"] <- "V2_LLUL1"
names(L_test)[names(L_test) == "V2"] <- "V2_LLUL2"


LL_time_size <- ggplot(L_test, aes(x = size, y = V2_LLtime)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "chromosome size", y = "Covariance", title="LL") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 10)
  )

LL_time_size

LL_UL_1_size <- ggplot(L_test, aes(x = size, y = V2_LLUL1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "chromosome size", y = "Covariance", title="LL_UL_1") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 10)
  )

LL_UL_1_size

LL_UL_2_size <- ggplot(L_test, aes(x = size, y = V2_LLUL2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "chromosome size", y = "Covariance", title="LL_UL_2") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 10))

LL_UL_2_size

png("Documents/FIBR_temporal/cvtk/chr_size_ex.png",  width=12, height=6, units = "in", res = 300)
plot_grid(CA_time_size, CA_TA_1_size, CA_TA_2_size, LL_time_size, LL_UL_1_size, LL_UL_2_size, nrow=2)
dev.off()
