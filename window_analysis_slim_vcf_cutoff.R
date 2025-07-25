library(ggplot2)
library(cowplot)
TA <- read.table("Documents/FIBR_temporal/slims_test/TA_vcf_slim_winds.txt", header=T)
head(TA)
summary(TA)
TA <- subset(TA, site_count < 1000)
q1 <- quantile(TA$mean_sel1, 0.95, na.rm = TRUE)
q2 <- quantile(TA$mean_sel2, 0.95, na.rm = TRUE)
qtot <- quantile(TA$mean_seltot, 0.95, na.rm = TRUE)
q1
q2
qtot

TA_obs <- read.table("Documents/FIBR_temporal/TA_99straight_cutoff_50kb.txt", header=T)
summary(TA_obs)
TA_obs <- subset(TA_obs, site_count < 1000)
TA_obs$ID <- paste0(TA_obs$chr, "_", TA_obs$BP1)
TA_outs_sel1 <- subset(TA_obs, mean_sel1 > q1)
nrow(TA_outs_sel1)
(nrow(TA_outs_sel1)/nrow(TA_obs))*100
TA_outs_sel2 <- subset(TA_obs, mean_sel2 > q2)
nrow(TA_outs_sel2)
(nrow(TA_outs_sel2)/nrow(TA_obs))*100
TA_outs_seltot <- subset(TA_obs, mean_seltot > qtot)
nrow(TA_outs_seltot)
(nrow(TA_outs_seltot)/nrow(TA_obs))*100


#CA
CA <- read.table("Documents/FIBR_temporal/slims_test/CA_vcf_slim_winds.txt", header=T)
summary(CA)
summary(TA)
CA <- subset(CA, site_count < 1000)
CAq1 <- quantile(CA$mean_sel1, 0.95, na.rm = TRUE)
CAq2 <- quantile(CA$mean_sel2, 0.95, na.rm = TRUE)
CAqtot <- quantile(CA$mean_seltot, 0.95, na.rm = TRUE)
CAq1
CAq2
CAqtot

CA_obs <- read.table("Documents/FIBR_temporal/CA_99straight_cutoff_50kb.txt", header=T)
CA_obs <- subset(CA_obs, site_count < 1000)
CA_obs$ID <- paste0(CA_obs$chr, "_", CA_obs$BP1)
CA_outs_sel1 <- subset(CA_obs, mean_sel1 > CAq1)
nrow(CA_outs_sel1)
(nrow(CA_outs_sel1)/nrow(CA_obs))*100
CA_outs_sel2 <- subset(CA_obs, mean_sel2 > CAq2)
nrow(CA_outs_sel2)
(nrow(CA_outs_sel2)/nrow(CA_obs))*100
CA_outs_seltot <- subset(CA_obs, mean_seltot > CAqtot)
nrow(CA_outs_seltot)
(nrow(CA_outs_seltot)/nrow(CA_obs))*100

#UL
UL <- read.table("Documents/FIBR_temporal/slims_test/UL_vcf_slim_winds.txt", header=T)
summary(UL)
summary(TA)
UL <- subset(UL, site_count < 1000)

ULq1 <- quantile(UL$mean_sel1, 0.95, na.rm = TRUE)
ULq2 <- quantile(UL$mean_sel2, 0.95, na.rm = TRUE)
ULqtot <- quantile(UL$mean_seltot, 0.95, na.rm = TRUE)
ULq1
ULq2
ULqtot

UL_obs <- read.table("Documents/FIBR_temporal/UL_99straight_cutoff_50kb.txt", header=T)
UL_obs <- subset(UL_obs, site_count < 1000)
UL_obs$ID <- paste0(UL_obs$chr, "_", UL_obs$BP1)
UL_outs_sel1 <- subset(UL_obs, mean_sel1 > ULq1)
nrow(UL_outs_sel1)
(nrow(UL_outs_sel1)/nrow(UL_obs))*100
UL_outs_sel2 <- subset(UL_obs, mean_sel2 > ULq2)
nrow(UL_outs_sel2)
(nrow(UL_outs_sel2)/nrow(UL_obs))*100
UL_outs_seltot <- subset(UL_obs, mean_seltot > ULqtot)
nrow(UL_outs_seltot)
(nrow(UL_outs_seltot)/nrow(UL_obs))*100

#LL
LL <- read.table("Documents/FIBR_temporal/slims_test/LL_vcf_slim_winds.txt", header=T)
summary(LL)
summary(TA)
LL <- subset(UL, site_count < 1000)

LLq1 <- quantile(LL$mean_sel1, 0.95, na.rm = TRUE)
LLq2 <- quantile(LL$mean_sel2, 0.95, na.rm = TRUE)
LLqtot <- quantile(LL$mean_seltot, 0.95, na.rm = TRUE)
LLq1
LLq2
LLqtot

LL_obs <- read.table("Documents/FIBR_temporal/LL_99straight_cutoff_50kb.txt", header=T)
LL_obs <- subset(LL_obs, site_count < 1000)
LL_obs$ID <- paste0(LL_obs$chr, "_", LL_obs$BP1)
LL_outs_sel1 <- subset(LL_obs, mean_sel1 > LLq1)
nrow(LL_outs_sel1)
(nrow(LL_outs_sel1)/nrow(LL_obs))*100
LL_outs_sel2 <- subset(LL_obs, mean_sel2 > LLq2)
nrow(LL_outs_sel2)
(nrow(LL_outs_sel2)/nrow(LL_obs))*100
LL_outs_seltot <- subset(LL_obs, mean_seltot > LLqtot)
nrow(LL_outs_seltot)
(nrow(LL_outs_seltot)/nrow(LL_obs))*100

##explore chr sel 2
LL_chr <- subset(LL_obs, chr == "chr20")

TA_chr <- subset(TA_obs, chr == "chr20")
head(TA_chr)
CA_chr <- subset(CA_obs, chr == "chr20")

UL_chr <- subset(UL_obs, chr == "chr20")
head(TA_chr)
highlight_region <- data.frame(
  xmin = 9750001,
  xmax = 9950001,
  ymin = -Inf,
  ymax = Inf
)

TA_chrom_plot <- ggplot(TA_chr, aes(x=(BP1+500), y=mean_sel2))+
  geom_point(aes(color = mean_sel2 > q2), size = 0.5)+
#  geom_rect(data = highlight_region,
#            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0,0.3))+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

TA_chrom_plot

LL_chrom_plot <- ggplot(LL_chr, aes(x=(BP1+500), y=mean_sel2))+
  geom_point(aes(color = mean_sel2 > LLq2), size = 0.5)+
#  geom_rect(data = highlight_region,
#            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

LL_chrom_plot

UL_chrom_plot <- ggplot(UL_chr, aes(x=(BP1+500), y=mean_sel2))+
  geom_point(aes(color = mean_sel2 > ULq2), size = 0.5)+
#  geom_rect(data = highlight_region,
#            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

UL_chrom_plot

CA_chrom_plot <- ggplot(CA_chr, aes(x=(BP1+500), y=mean_sel2))+
  geom_point(aes(color = mean_sel2 > CAq2), size = 0.5)+
#  geom_rect(data = highlight_region,
#            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

CA_chrom_plot

plot_grid(CA_chrom_plot, TA_chrom_plot, UL_chrom_plot, LL_chrom_plot)


##explore chr sel 1
LL_chr <- subset(LL_obs, chr == "chr15")
TA_chr <- subset(TA_obs, chr == "chr15")
CA_chr <- subset(CA_obs, chr == "chr15")
UL_chr <- subset(UL_obs, chr == "chr15")

head(TA_chr)

highlight_region <- data.frame(
  xmin = 9750001,
  xmax = 9950001,
  ymin = -Inf,
  ymax = Inf
)
head(TA_chr)

TA_chrom_plot <- ggplot(TA_chr, aes(x=(BP1+500), y=mean_sel1))+
  geom_point(aes(color = mean_sel1 > q1), size = 0.5)+
  #  geom_rect(data = highlight_region,
  #            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0,0.3))+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

TA_chrom_plot

LL_chrom_plot <- ggplot(LL_chr, aes(x=(BP1+500), y=mean_sel1))+
  geom_point(aes(color = mean_sel1 > LLq1), size = 0.5)+
  #  geom_rect(data = highlight_region,
  #            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

LL_chrom_plot

UL_chrom_plot <- ggplot(UL_chr, aes(x=(BP1+500), y=mean_sel1))+
  geom_point(aes(color = mean_sel1 > ULq2), size = 0.5)+
  #  geom_rect(data = highlight_region,
  #            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

UL_chrom_plot

CA_chrom_plot <- ggplot(CA_chr, aes(x=(BP1+500), y=mean_sel1))+
  geom_point(aes(color = mean_sel1 > CAq2), size = 0.5)+
  #  geom_rect(data = highlight_region,
  #            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #            fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+theme(legend.position = "none")

CA_chrom_plot

plot_grid(CA_chrom_plot, TA_chrom_plot, UL_chrom_plot, LL_chrom_plot)

##sup figure
head(TA)
head(TA_obs)

TA_tmp <- data.frame(mean_sel1 = TA$mean_sel1, dataset= "simulated")
head(TA_tmp)
TA_obs_tmp <- data.frame(mean_sel1 = TA_obs$mean_sel1, dataset= "observed")
final_df <- data.frame(rbind(TA_tmp, TA_obs_tmp))
head(final_df)                       

sup_fig <- ggplot(final_df, aes(x = dataset, y = mean_sel1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  labs(
    x = "Dataset",
    y = "Change in Allele Frequency (ΔAF)",
    fill = "Dataset",
    title = "Comparison of ΔAF Distributions by Starting Frequency Bin"
  ) +
  theme_minimal()


sup_fig


##overlap analysis
# Combine all IDs
all_sel1_out_ID <- unique(c(CA_outs_sel1$ID, TA_outs_sel1$ID, UL_outs_sel1$ID, LL_outs_sel1$ID))


# Create a data frame with all IDs
id_df <- data.frame(ID = all_sel1_out_ID)
head(id_df)

# Count presence in each dataframe
id_df$in_CA <- id_df$ID %in% CA_outs_sel1$ID
id_df$in_TA <- id_df$ID %in% TA_outs_sel1$ID
id_df$in_UL <- id_df$ID %in% UL_outs_sel1$ID
id_df$in_LL <- id_df$ID %in% LL_outs_sel1$ID

# Sum across rows to count how many dataframes each ID is in
id_df$count <- rowSums(id_df[ , c("in_CA", "in_TA", "in_UL", "in_LL")])
head(id_df)
# Now you can filter:
in_all_4 <- subset(id_df, count == 4)
nrow(in_all_4)
in_all_4
in_any_3 <- subset(id_df, count == 3)
nrow(in_any_3)
in_any_3
in_any_2 <- subset(id_df, count == 2)
nrow(in_any_2)
in_only_1 <- subset(id_df, count == 1)
nrow(in_only_1)

##in sel2

# Combine all IDs
all_sel2_out_ID <- unique(c(CA_outs_sel2$ID, TA_outs_sel2$ID, UL_outs_sel2$ID, LL_outs_sel2$ID))


# Create a data frame with all IDs
id_df_2 <- data.frame(ID = all_sel2_out_ID)

# Count presence in each dataframe
id_df_2$in_CA <- id_df_2$ID %in% CA_outs_sel2$ID
id_df_2$in_TA <- id_df_2$ID %in% TA_outs_sel2$ID
id_df_2$in_UL <- id_df_2$ID %in% UL_outs_sel2$ID
id_df_2$in_LL <- id_df_2$ID %in% LL_outs_sel2$ID

# Sum across rows to count how many dataframes each ID is in
id_df_2$count <- rowSums(id_df_2[ , c("in_CA", "in_TA", "in_UL", "in_LL")])
head(id_df_2)
# Now you can filter:
in_all_4 <- subset(id_df_2, count == 4)
nrow(in_all_4)
in_all_4
in_any_3 <- subset(id_df_2, count == 3)
nrow(in_any_3)
in_any_2 <- subset(id_df_2, count == 2)
nrow(in_any_2)
in_only_1 <- subset(id_df_2, count == 1)
nrow(in_only_1)

head(CA)

# Combine all IDs
all_seltot_out_ID <- unique(c(CA_outs_seltot$ID, TA_outs_seltot$ID, UL_outs_seltot$ID, LL_outs_seltot$ID))


# Create a data frame with all IDs
id_df_tot <- data.frame(ID = all_seltot_out_ID)

# Count presence in each dataframe
id_df_tot$in_CA <- id_df_tot$ID %in% CA_outs_seltot$ID
id_df_tot$in_TA <- id_df_tot$ID %in% TA_outs_seltot$ID
id_df_tot$in_UL <- id_df_tot$ID %in% UL_outs_seltot$ID
id_df_tot$in_LL <- id_df_tot$ID %in% LL_outs_seltot$ID

# Sum across rows to count how many dataframes each ID is in
id_df_tot$count <- rowSums(id_df_tot[ , c("in_CA", "in_TA", "in_UL", "in_LL")])
head(id_df_tot)
# Now you can filter:
in_all_4 <- subset(id_df_tot, count == 4)
nrow(in_all_4)
in_all_4
in_any_3 <- subset(id_df_tot, count == 3)
nrow(in_any_3)
in_any_3
in_any_2 <- subset(id_df_tot, count == 2)
nrow(in_any_2)
in_only_1 <- subset(id_df_tot, count == 1)
nrow(in_only_1)


