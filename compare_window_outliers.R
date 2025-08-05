TA <- read.table("Documents/FIBR_temporal/TA_99straight_cutoff_50kb.txt", header=T)
summary(TA$site_count)
TA <- subset(TA, site_count < 1000)
CA <- read.table("Documents/FIBR_temporal/CA_99straight_cutoff_50kb.txt", header=T)
CA <- subset(CA, site_count < 1000)
UL <- read.table("Documents/FIBR_temporal/UL_99straight_cutoff_50kb.txt", header=T)
UL <- subset(UL, site_count < 1000)
LL <- read.table("Documents/FIBR_temporal/LL_99straight_cutoff_50kb.txt", header=T)
LL <- subset(LL, site_count < 1000)

head(CA)
CA$ID <- paste0(CA$chr, "_", CA$BP1)
TA$ID <- paste0(TA$chr, "_", TA$BP1)
UL$ID <- paste0(UL$chr, "_", UL$BP1)
LL$ID <- paste0(LL$chr, "_", LL$BP1)

summary(CA$outlier_count_sel1/CA$site_count)

##in sel1
#CA_sel1_out <- subset(CA, outlier_count_sel1/site_count > 100)
CA_sel1_out <- subset(CA, (CA$outlier_count_sel1/CA$site_count) > 0.2)
nrow(CA_sel1_out)
TA_sel1_out <- subset(TA, (TA$outlier_count_sel1/TA$site_count) > 0.2)
nrow(TA_sel1_out)
UL_sel1_out <- subset(UL,  (UL$outlier_count_sel1/UL$site_count) > 0.2)
nrow(UL_sel1_out)
LL_sel1_out <- subset(LL, (LL$outlier_count_sel1/LL$site_count) > 0.2)
nrow(LL_sel1_out)

# Combine all IDs
all_sel1_out_ID <- unique(c(CA_sel1_out$ID, TA_sel1_out$ID, UL_sel1_out$ID, LL_sel1_out$ID))


# Create a data frame with all IDs
id_df <- data.frame(ID = all_sel1_out_ID)
head(id_df)

# Count presence in each dataframe
id_df$in_CA <- id_df$ID %in% CA_sel1_out$ID
id_df$in_TA <- id_df$ID %in% TA_sel1_out$ID
id_df$in_UL <- id_df$ID %in% UL_sel1_out$ID
id_df$in_LL <- id_df$ID %in% LL_sel1_out$ID

# Sum across rows to count how many dataframes each ID is in
id_df$count <- rowSums(id_df[ , c("in_CA", "in_TA", "in_UL", "in_LL")])
head(id_df)
# Now you can filter:
in_all_4 <- subset(id_df, count == 4)
nrow(in_all_4)
in_all_4
in_any_3 <- subset(id_df, count == 3)
nrow(in_any_3)
in_any_2 <- subset(id_df, count == 2)
nrow(in_any_2)
in_only_1 <- subset(id_df, count == 1)
nrow(in_only_1)

head(CA)

##in sel2
CA_sel2_out <- subset(CA,  (CA$outlier_count_sel2/CA$site_count) > 0.2)
nrow(CA_sel2_out)
TA_sel2_out <- subset(TA, (TA$outlier_count_sel2/TA$site_count) > 0.2)
nrow(TA_sel2_out)
UL_sel2_out <- subset(UL, (UL$outlier_count_sel2/UL$site_count) > 0.2)
nrow(UL_sel2_out)
LL_sel2_out <- subset(LL, (LL$outlier_count_sel2/LL$site_count) > 0.2)
nrow(LL_sel2_out)

# Combine all IDs
all_sel2_out_ID <- unique(c(CA_sel2_out$ID, TA_sel2_out$ID, UL_sel2_out$ID, LL_sel2_out$ID))


# Create a data frame with all IDs
id_df_2 <- data.frame(ID = all_sel2_out_ID)

# Count presence in each dataframe
id_df_2$in_CA <- id_df_2$ID %in% CA_sel2_out$ID
id_df_2$in_TA <- id_df_2$ID %in% TA_sel2_out$ID
id_df_2$in_UL <- id_df_2$ID %in% UL_sel2_out$ID
id_df_2$in_LL <- id_df_2$ID %in% LL_sel2_out$ID

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
##in sel_tot
CA_seltot_out <- subset(CA,  (CA$outlier_count_seltot/CA$site_count) > 0.2)
nrow(CA_seltot_out)
TA_seltot_out <- subset(TA, (TA$outlier_count_seltot/TA$site_count) > 0.2)
nrow(TA_seltot_out)
UL_seltot_out <- subset(UL, (UL$outlier_count_seltot/UL$site_count) > 0.2)
nrow(UL_seltot_out)
LL_seltot_out <- subset(LL, (LL$outlier_count_seltot/LL$site_count) > 0.2)
nrow(LL_seltot_out)

# Combine all IDs
all_seltot_out_ID <- unique(c(CA_seltot_out$ID, TA_seltot_out$ID, UL_seltot_out$ID, LL_seltot_out$ID))


# Create a data frame with all IDs
id_df_tot <- data.frame(ID = all_seltot_out_ID)

# Count presence in each dataframe
id_df_tot$in_CA <- id_df_tot$ID %in% CA_seltot_out$ID
id_df_tot$in_TA <- id_df_tot$ID %in% TA_seltot_out$ID
id_df_tot$in_UL <- id_df_tot$ID %in% UL_seltot_out$ID
id_df_tot$in_LL <- id_df_tot$ID %in% LL_seltot_out$ID

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

####explore
CA_chr <- subset(CA, chr == "chr15")
TA_chr <- subset(TA, chr == "chr15")
UL_chr <- subset(UL, chr == "chr15")
LL_chr <- subset(LL, chr == "chr15")

head(CA_chr)

CA_chrom_plot <- ggplot(CA_chr, aes(x=(BP1+500), y=(outlier_count_sel1/site_count)))+
  geom_point(size = 0.5)+
#  geom_point(aes(x = BP1 + 500, y = outlier_count_sel2), color = "red", size = 0.5)+ 
#  geom_point(aes(x = BP1 + 500, y = outlier_count_seltot), color = "blue", size = 0.5)+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()

CA_chrom_plot

TA_chrom_plot <- ggplot(TA_chr, aes(x=(BP1+500), y=(outlier_count_sel1/site_count)))+
  geom_point(size = 0.5)+
#  geom_point(aes(x = BP1 + 500, y = outlier_count_sel2), color = "red", size = 0.5)+ 
#  geom_point(aes(x = BP1 + 500, y = outlier_count_seltot), color = "blue", size = 0.5)+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()

TA_chrom_plot

UL_chrom_plot <- ggplot(UL_chr, aes(x=(BP1+500), y=(outlier_count_sel1/site_count)))+
  geom_point(size = 0.5)+
 # geom_point(aes(x = BP1 + 500, y = outlier_count_sel2), color = "red", size = 0.5)+ 
#  geom_point(aes(x = BP1 + 500, y = outlier_count_seltot), color = "blue", size = 0.5)+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()

UL_chrom_plot

LL_chrom_plot <- ggplot(LL_chr, aes(x=(BP1+500), y=(outlier_count_sel1/site_count)))+
  geom_point(size = 0.5)+
#  geom_point(aes(x = BP1 + 500, y = outlier_count_sel2), color = "red", size = 0.5)+ 
#  geom_point(aes(x = BP1 + 500, y = outlier_count_seltot), color = "blue", size = 0.5)+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()

LL_chrom_plot

plot_grid(CA_chrom_plot, TA_chrom_plot, LL_chrom_plot, UL_chrom_plot)


##make genome-wide plots for each measure for sup###
