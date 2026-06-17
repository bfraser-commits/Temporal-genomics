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

###sup figures
head(TA)
TA_obs_tmp <- data.frame(mean_sel1 = TA_obs$mean_sel1, mean_sel2 = TA_obs$mean_sel2, mean_seltot = TA_obs$mean_seltot, dataset= "observed")
TA_slim_tmp <-  data.frame(mean_sel1 = TA$mean_sel1, mean_sel2 = TA$mean_sel2, mean_seltot = TA$mean_seltot, dataset= "simulated")

TA_final_df <- data.frame(rbind(TA_obs_tmp, TA_slim_tmp))
head(final_df)


sel_1_fig <- ggplot(TA_final_df, aes(x = dataset, y = mean_sel1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = q1, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 1",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

sel_1_fig

sel_2_fig <- ggplot(TA_final_df, aes(x = dataset, y = mean_sel2, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = q2, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 2",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

sel_2_fig


sel_tot_fig <- ggplot(TA_final_df, aes(x = dataset, y = mean_seltot, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = qtot, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Total Period",
    x = "Dataset",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

sel_tot_fig

combined_plot <- plot_grid(sel_1_fig, sel_2_fig, sel_tot_fig, ncol=1, labels = c("A", "B", "C"), align = "v")

final_plot <- plot_grid(
  ggdraw() + 
    draw_label("TA",  fontface = 'bold', size = 16, x = 0.5, hjust = 0.5), 
  combined_plot, ncol = 1,  rel_heights = c(0.1, 1.5)  
)

final_plot

pdf("Documents/FIBR_temporal/slims_test/TA_sup_fig_window.pdf", height = 6, width = 3)
final_plot
dev.off()

head(CA)
CA_obs_tmp <- data.frame(mean_sel1 = CA_obs$mean_sel1, mean_sel2 = CA_obs$mean_sel2, mean_seltot = CA_obs$mean_seltot, dataset= "observed")
CA_slim_tmp <-  data.frame(mean_sel1 = CA$mean_sel1, mean_sel2 = CA$mean_sel2, mean_seltot = CA$mean_seltot, dataset= "simulated")

CA_final_df <- data.frame(rbind(CA_obs_tmp, CA_slim_tmp))
head(final_df)


CA_sel_1_fig <- ggplot(CA_final_df, aes(x = dataset, y = mean_sel1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = CAq1, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 1",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

CA_sel_1_fig

CA_sel_2_fig <- ggplot(CA_final_df, aes(x = dataset, y = mean_sel2, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = CAq2, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 2",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

CA_sel_2_fig


CA_sel_tot_fig <- ggplot(CA_final_df, aes(x = dataset, y = mean_seltot, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = CAqtot, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Total Period",
    x = "Dataset",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

CA_sel_tot_fig

CA_combined_plot <- plot_grid(CA_sel_1_fig, CA_sel_2_fig, CA_sel_tot_fig, ncol=1, labels = c("A", "B", "C"), align = "v")

CA_final_plot <- plot_grid(
  ggdraw() + 
    draw_label("CA",  fontface = 'bold', size = 16, x = 0.5, hjust = 0.5), 
  CA_combined_plot, ncol = 1,  rel_heights = c(0.1, 1.5)  
)



pdf("Documents/FIBR_temporal/slims_test/CA_sup_fig_window.pdf", height = 6, width = 3)
CA_final_plot
dev.off()


head(UL)
UL_obs_tmp <- data.frame(mean_sel1 = UL_obs$mean_sel1, mean_sel2 = UL_obs$mean_sel2, mean_seltot = UL_obs$mean_seltot, dataset= "observed")
UL_slim_tmp <-  data.frame(mean_sel1 = UL$mean_sel1, mean_sel2 = UL$mean_sel2, mean_seltot = UL$mean_seltot, dataset= "simulated")

UL_final_df <- data.frame(rbind(UL_obs_tmp, UL_slim_tmp))
head(UL_final_df)


UL_sel_1_fig <- ggplot(UL_final_df, aes(x = dataset, y = mean_sel1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = ULq1, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 1",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

UL_sel_1_fig
ULq2

UL_sel_2_fig <- ggplot(UL_final_df, aes(x = dataset, y = mean_sel2, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = ULq2, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 2",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

UL_sel_2_fig

ULqtot

UL_sel_tot_fig <- ggplot(UL_final_df, aes(x = dataset, y = mean_seltot, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = ULqtot, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Total Period",
    x = "Dataset",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

UL_sel_tot_fig

UL_combined_plot <- plot_grid(UL_sel_1_fig, UL_sel_2_fig, UL_sel_tot_fig, ncol=1, labels = c("A", "B", "C"), align = "v")

UL_final_plot <- plot_grid(
  ggdraw() + 
    draw_label("UL",  fontface = 'bold', size = 16, x = 0.5, hjust = 0.5), 
  UL_combined_plot, ncol = 1,  rel_heights = c(0.1, 1.5)  
)



pdf("Documents/FIBR_temporal/slims_test/UL_sup_fig_window.pdf", height = 6, width = 3)
UL_final_plot
dev.off()


head(LL)
LL_obs_tmp <- data.frame(mean_sel1 = LL_obs$mean_sel1, mean_sel2 = LL_obs$mean_sel2, mean_seltot = LL_obs$mean_seltot, dataset= "observed")
LL_slim_tmp <-  data.frame(mean_sel1 = LL$mean_sel1, mean_sel2 = LL$mean_sel2, mean_seltot = LL$mean_seltot, dataset= "simulated")

LL_final_df <- data.frame(rbind(LL_obs_tmp, LL_slim_tmp))
head(LL_final_df)

LLq1

LL_sel_1_fig <- ggplot(LL_final_df, aes(x = dataset, y = mean_sel1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = LLq1, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 1",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

LL_sel_1_fig
LLq2

LL_sel_2_fig <- ggplot(LL_final_df, aes(x = dataset, y = mean_sel2, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = LLq2, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Period 2",
    x = "",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

LL_sel_2_fig

LLqtot

LL_sel_tot_fig <- ggplot(LL_final_df, aes(x = dataset, y = mean_seltot, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  geom_hline(yintercept = LLqtot, linetype = "dashed", color = "darkgrey")+
  labs(
    title = "Total Period",
    x = "Dataset",
    y = "Mean delta AF",
    fill = "Dataset",
  ) +
  theme_minimal()+
  theme(legend.position = "none")

LL_sel_tot_fig

LL_combined_plot <- plot_grid(LL_sel_1_fig, LL_sel_2_fig, LL_sel_tot_fig, ncol=1, labels = c("A", "B", "C"), align = "v")

LL_final_plot <- plot_grid(
  ggdraw() + 
    draw_label("LL",  fontface = 'bold', size = 16, x = 0.5, hjust = 0.5), 
  LL_combined_plot, ncol = 1,  rel_heights = c(0.1, 1.5)  
)



pdf("Documents/FIBR_temporal/slims_test/LL_sup_fig_window.pdf", height = 6, width = 3)
LL_final_plot
dev.off()

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

candidate_sel1 <- subset(id_df, count >= 3)
head(candidate_sel1)
nrow(candidate_sel1)
write.table(candidate_sel1, "Documents/FIBR_temporal/slims_test/candidate_window_sel1.txt", sep = "\t", row.names=F)
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
in_any_3
in_any_2 <- subset(id_df_2, count == 2)
nrow(in_any_2)
in_only_1 <- subset(id_df_2, count == 1)
nrow(in_only_1)

head(CA)
candidate_sel2 <- subset(id_df_2, count >=3)
nrow(candidate_sel2)
write.table(candidate_sel2, "Documents/FIBR_temporal/slims_test/candidate_window_sel2.txt", sep = "\t", row.names=F)

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

candidate_regions_total <- subset(id_df_tot, count >= 3)
nrow(candidate_regions_total)
write.table(candidate_regions_total, file = "Documents/FIBR_temporal/slims_test/candidate_window_total.txt", sep = "\t", quote=F, row.names=F)

###make upset plots
library(UpSetR)
head(id_df)
up_id_df <- data.frame("CA" = id_df$in_CA, "TA" = id_df$in_TA, "LL" = id_df$in_LL, "UL" = id_df$in_UL)
up_id_df <- as.data.frame(lapply(up_id_df, as.numeric))
head(up_id_df)

# Create UpSet plot
sel_1_upset <- upset(up_id_df,
      sets = c("CA", "TA", "LL", "UL"),
      #order.by = "freq",
      keep.order = TRUE, 
      sets.bar.color = c("CA" = "#FFB90F",
                         "TA" = "#FF7F00",
                         "LL" = "#5D478B",
                         "UL" = "#CD0000"))
    


sel_1_upset

head(id_df_2)
up_id_df_2 <- data.frame("CA" = id_df_2$in_CA, "TA" = id_df_2$in_TA, "LL" = id_df_2$in_LL, "UL" = id_df_2$in_UL)
up_id_df_2 <- as.data.frame(lapply(up_id_df_2, as.numeric))
head(up_id_df_2)

sel_2_upset <- upset(up_id_df_2,
                     sets = c("CA", "TA", "LL", "UL"),
                     #order.by = "freq",
                     keep.order = TRUE,
                     sets.bar.color = c("CA" = "#FFB90F",
                                        "TA" = "#FF7F00",
                                        "LL" = "#5D478B",
                                        "UL" = "#CD0000"))
sel_2_upset

head(id_df_tot)
up_id_df_tot <- data.frame("CA" = id_df_tot$in_CA, "TA" = id_df_tot$in_TA, "LL" = id_df_tot$in_LL, "UL" = id_df_tot$in_UL)
up_id_df_tot <- as.data.frame(lapply(up_id_df_tot, as.numeric))
head(up_id_df_tot)

sel_tot_upset <- upset(up_id_df_tot,
                     sets = c("CA", "TA", "LL", "UL"),
                     #order.by = "freq",
                     keep.order = TRUE,
                     sets.bar.color = c("CA" = "#FFB90F",
                                        "TA" = "#FF7F00",
                                        "LL" = "#5D478B",
                                        "UL" = "#CD0000"),
                     title = "YOUR TITLE HERE")
sel_tot_upset

pdf("Documents/FIBR_temporal/slims_test/sel_1_upset.pdf", width=5, height=5)
sel_1_upset
dev.off()
pdf("Documents/FIBR_temporal/slims_test/sel_2_upset.pdf", width=5, height=5)
sel_2_upset
dev.off()
pdf("Documents/FIBR_temporal/slims_test/sel_tot_upset.pdf", width=5, height=5)
sel_tot_upset
dev.off()
