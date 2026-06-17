library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(cowplot)

#take candidate loci for each time period
sel1 <- read.table("Documents/FIBR_temporal/slims_test/candidate_window_sel1.txt", header=T)
split_vals <- do.call(rbind, strsplit(sel1$ID, "_"))
sel1$chr <- split_vals[, 1]
sel1$pos <- as.numeric(split_vals[, 2])
sel1$end <- sel1$pos + 49999
head(sel1)

sel2 <- read.table("Documents/FIBR_temporal/slims_test/candidate_window_sel2.txt", header=T)
split_vals <- do.call(rbind, strsplit(sel2$ID, "_"))
sel2$chr <- split_vals[, 1]
sel2$pos <- as.numeric(split_vals[, 2])
sel2$end <- sel2$pos + 49999

seltot <- read.table("Documents/FIBR_temporal/slims_test/candidate_window_total.txt", header=T)
split_vals <- do.call(rbind, strsplit(seltot$ID, "_"))
seltot$chr <- split_vals[, 1]
seltot$pos <- as.numeric(split_vals[, 2])
seltot$end <- seltot$pos + 49999
head(seltot)

check_in_sel <- function(chr, pos, sel_df) {
  any(sel_df$chr == chr & sel_df$pos <= pos & sel_df$end >= pos)
}

#take AF for each pop

GH_AF <- data.frame(fread("Documents/FIBR_temporal/slims_test/GH.frq.txt"))
colnames(GH_AF) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
GH_AF$mutation_id <- paste0(GH_AF$chr,"_",GH_AF$pos)
GH_AF <- GH_AF[grep("chr",GH_AF$chr),]

CA_intro_AF_13 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/CA_2013.frq.txt")))
colnames(CA_intro_AF_13) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
CA_intro_AF_13$mutation_id <- paste0(CA_intro_AF_13$chr,"_",CA_intro_AF_13$pos)
CA_intro_AF_13 <- CA_intro_AF_13[grep("chr",CA_intro_AF_13$chr),]

CA_intro_AF_18 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/CA_2018.frq.txt")))
colnames(CA_intro_AF_18) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
CA_intro_AF_18$mutation_id <- paste0(CA_intro_AF_18$chr,"_",CA_intro_AF_18$pos)
CA_intro_AF_18 <- CA_intro_AF_18[grep("chr",CA_intro_AF_18$chr),]

#head(GH_AF)
CA_intro_tmp <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, p1 = GH_AF$REF_FREQ,
                           pmid =  CA_intro_AF_13$REF_FREQ, pt = CA_intro_AF_18$REF_FREQ)

head(CA_intro_tmp)

#switch to minor allele
#CA_intro_tmp$p1_new <- ifelse(CA_intro_tmp$p1 > 0.5, 1 - CA_intro_tmp$p1, CA_intro_tmp$p1)
#CA_intro_tmp$pmid_new <- ifelse(CA_intro_tmp$p1 > 0.5, 1 - CA_intro_tmp$pmid, CA_intro_tmp$pmid)
#CA_intro_tmp$pt_new <- ifelse(CA_intro_tmp$p1 > 0.5, 1 - CA_intro_tmp$pt, CA_intro_tmp$pt)

#get allele frequency change
CA_intro_tmp$sel_coef_1 <- (CA_intro_tmp$pmid - CA_intro_tmp$p1)
CA_intro_tmp$sel_coef_2 <- (CA_intro_tmp$pt - CA_intro_tmp$pmid)
CA_intro_tmp$sel_coef_total <- (CA_intro_tmp$pt - CA_intro_tmp$p1)

#get only windows in that pop
CA_sel1 <- subset(sel1, in_CA == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
CA_intro_tmp$in_sel1 <- mapply(check_in_sel, CA_intro_tmp$chr, CA_intro_tmp$pos, MoreArgs = list(sel_df = CA_sel1))

#then just subset the alleles that are under selection
#update to include cut-offs from slims
CA_intro_tmp$in_sel1_allele <- ifelse(CA_intro_tmp$sel_coef_1 > 0.21 | CA_intro_tmp$sel_coef_1 < -0.20, TRUE, FALSE)
CA_intro_tmp$both_in_sel1 <- CA_intro_tmp$in_sel1 & CA_intro_tmp$in_sel1_allele

CA_sel_1_test <- subset(CA_intro_tmp, both_in_sel1 == "TRUE")

head(CA_sel_1_test)
nrow(CA_sel_1_test)
summary(CA_sel_1_test$sel_coef_1)

CA_sel1_plot <- ggplot(CA_sel_1_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in 1st Time Interval",
    x = "Delta AF Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)

CA_sel1_plot

summary(CA_sel_1_test)
CA_sel_1_test_pos <- subset(CA_sel_1_test, sel_coef_1 > 0)
quantile(CA_sel_1_test_pos$sel_coef_2, probs = c(0.05, 0.5, 0.95))
CA_sel_1_test_neg <- subset(CA_sel_1_test, sel_coef_1 < 0)
quantile(CA_sel_1_test_neg$sel_coef_2, probs = c(0.05, 0.5, 0.95))


#Now do sel2
#get only windows in that pop
CA_sel2 <- subset(sel2, in_CA == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
CA_intro_tmp$in_sel2 <- mapply(check_in_sel, CA_intro_tmp$chr, CA_intro_tmp$pos, MoreArgs = list(sel_df = CA_sel2))
head(CA_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
CA_intro_tmp$in_sel2_allele <- ifelse(CA_intro_tmp$sel_coef_2 > 0.20 | CA_intro_tmp$sel_coef_2 < - 0.20, TRUE, FALSE)
CA_intro_tmp$both_in_sel2 <- CA_intro_tmp$in_sel2 & CA_intro_tmp$in_sel2_allele

CA_sel_2_test <- subset(CA_intro_tmp, both_in_sel2 == "TRUE")

head(CA_sel_1_test)
nrow(CA_sel_2_test)

CA_sel2_plot <- ggplot(CA_sel_2_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in 2nd Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)


CA_sel2_plot

CA_sel_2_test_pos <- subset(CA_sel_2_test, sel_coef_2 > 0)
quantile(CA_sel_2_test_pos$sel_coef_1, probs = c(0.05, 0.5, 0.95))
CA_sel_2_test_neg <- subset(CA_sel_2_test, sel_coef_2 < 0)
quantile(CA_sel_2_test_neg$sel_coef_1, probs = c(0.05, 0.5, 0.95))


#now total plot
#get only windows in that pop
CA_seltot <- subset(seltot, in_CA == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
CA_intro_tmp$in_seltot <- mapply(check_in_sel, CA_intro_tmp$chr, CA_intro_tmp$pos, MoreArgs = list(sel_df = CA_seltot))
head(CA_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
CA_intro_tmp$in_seltot_allele <- ifelse(CA_intro_tmp$sel_coef_total > 0.24 | CA_intro_tmp$sel_coef_total < - 0.23, TRUE, FALSE)
CA_intro_tmp$both_in_seltot <- CA_intro_tmp$in_seltot & CA_intro_tmp$in_seltot_allele

CA_sel_tot_test <- subset(CA_intro_tmp, both_in_seltot == "TRUE")

head(CA_sel_1_test)
nrow(CA_sel_tot_test)

CA_seltot_plot <- ggplot(CA_sel_tot_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in Total Time",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)


CA_seltot_plot

CA_sel_tot_pos_1 <- subset(CA_sel_tot_test, sel_coef_total > 0)
quantile(CA_sel_2_test_pos$sel_coef_1, probs = c(0.05, 0.5, 0.95))
quantile(CA_sel_2_test_pos$sel_coef_2, probs = c(0.05, 0.5, 0.95))

CA_sel_tot_neg_1 <- subset(CA_sel_tot_test, sel_coef_total < 0)
quantile(CA_sel_2_test_neg$sel_coef_1, probs = c(0.05, 0.5, 0.95))




#CA_combined <- plot_grid(CA_sel1_plot, CA_sel2_plot, CA_seltot_plot, ncol=3, labels = c('(A)', '(B)', '(C)'))
CA_combined <- plot_grid(CA_sel1_plot, CA_sel2_plot, CA_seltot_plot, ncol=3, labels = c('(D)', '(E)', '(F)'))

CA_final_plot <- ggdraw() +
  draw_label("CA", fontface = 'bold', size = 16, x = 0.5, y = 0.97, hjust = 0.5) +
  draw_plot(CA_combined, y = 0, height = 0.95)

png("Documents/FIBR_temporal/allele_trajectory/CA_AF_density_plot.png", width = 10, height = 4, units = "in", res = 300)
CA_final_plot
dev.off()




#######now do it for Taylor
TA_intro_AF_13 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/TA_2013.frq.txt")))
colnames(TA_intro_AF_13) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
TA_intro_AF_13$mutation_id <- paste0(TA_intro_AF_13$chr,"_",TA_intro_AF_13$pos)
TA_intro_AF_13 <- TA_intro_AF_13[grep("chr",TA_intro_AF_13$chr),]

TA_intro_AF_18 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/TA_2018.frq.txt")))
colnames(TA_intro_AF_18) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
TA_intro_AF_18$mutation_id <- paste0(TA_intro_AF_18$chr,"_",TA_intro_AF_18$pos)
TA_intro_AF_18 <- TA_intro_AF_18[grep("chr",TA_intro_AF_18$chr),]

#head(GH_AF)
TA_intro_tmp <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, p1 = GH_AF$REF_FREQ,
                           pmid =  TA_intro_AF_13$REF_FREQ, pt = TA_intro_AF_18$REF_FREQ)



#get allele frequency change
TA_intro_tmp$sel_coef_1 <- (TA_intro_tmp$pmid - TA_intro_tmp$p1)
TA_intro_tmp$sel_coef_2 <- (TA_intro_tmp$pt - TA_intro_tmp$pmid)
TA_intro_tmp$sel_coef_total <- (TA_intro_tmp$pt - TA_intro_tmp$p1)

head(sel1)
#get only windows in that pop
TA_sel1 <- subset(sel1, in_TA == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
TA_intro_tmp$in_sel1 <- mapply(check_in_sel, TA_intro_tmp$chr, TA_intro_tmp$pos, MoreArgs = list(sel_df = TA_sel1))

#then just subset the alleles that are under selection
#update to include cut-offs from slims
TA_intro_tmp$in_sel1_allele <- ifelse(TA_intro_tmp$sel_coef_1 > 0.21 | TA_intro_tmp$sel_coef_1 < -0.20, TRUE, FALSE)
TA_intro_tmp$both_in_sel1 <- TA_intro_tmp$in_sel1 & TA_intro_tmp$in_sel1_allele

TA_sel_1_test <- subset(TA_intro_tmp, both_in_sel1 == "TRUE")

#nrow(CA_sel_1_test)
nrow(TA_sel_1_test)

TA_sel1_plot <- ggplot(TA_sel_1_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  #  scale_fill_viridis_d(option = "C", direction = -1, na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  scale_fill_viridis_d(option = "C", na.value = "white") +
  theme_minimal()+
  labs(
    title = "Candidates in 1st Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)


TA_sel1_plot

#Now do sel2
#get only windows in that pop
TA_sel2 <- subset(sel2, in_TA == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
TA_intro_tmp$in_sel2 <- mapply(check_in_sel, TA_intro_tmp$chr, TA_intro_tmp$pos, MoreArgs = list(sel_df = TA_sel2))
head(TA_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
TA_intro_tmp$in_sel2_allele <- ifelse(TA_intro_tmp$sel_coef_2 > 0.20 | TA_intro_tmp$sel_coef_2 < -0.20, TRUE, FALSE)
TA_intro_tmp$both_in_sel2 <- TA_intro_tmp$in_sel2 & TA_intro_tmp$in_sel2_allele

TA_sel_2_test <- subset(TA_intro_tmp, both_in_sel2 == "TRUE")

nrow(CA_sel_2_test)
nrow(TA_sel_2_test)

TA_sel2_plot <- ggplot(TA_sel_2_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  #  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in 2nd Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)



TA_sel2_plot

#now total plot
#get only windows in that pop
TA_seltot <- subset(seltot, in_TA == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
TA_intro_tmp$in_seltot <- mapply(check_in_sel, TA_intro_tmp$chr, TA_intro_tmp$pos, MoreArgs = list(sel_df = TA_seltot))
head(TA_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
TA_intro_tmp$in_seltot_allele <- ifelse(TA_intro_tmp$sel_coef_total > 0.23 | TA_intro_tmp$sel_coef_total < - 0.22, TRUE, FALSE)
TA_intro_tmp$both_in_seltot <- TA_intro_tmp$in_seltot & TA_intro_tmp$in_seltot_allele

TA_sel_tot_test <- subset(TA_intro_tmp, both_in_seltot == "TRUE")

#nrow(CA_sel_tot_test)
nrow(TA_sel_tot_test)

TA_seltot_plot <- ggplot(TA_sel_tot_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  #  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in Total Time",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)



TA_seltot_plot

TA_combined_plot <- plot_grid(TA_sel1_plot, TA_sel2_plot, TA_seltot_plot, ncol=3, labels = c('(A)', '(B)', '(C)'))

TA_final_plot <- ggdraw() +
  draw_label("TA", fontface = 'bold', size = 16, x = 0.5, y = 0.97, hjust = 0.5) +
  draw_plot(TA_combined_plot, y = 0, height = 0.95)

png("Documents/FIBR_temporal/allele_trajectory/TA_AF_density_plot.png", width = 10, height = 4, units = "in", res = 300)
TA_final_plot
dev.off()

#######now do it for Upper Laj
UL_intro_AF_13 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/UL_2013.frq.txt")))
colnames(UL_intro_AF_13) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
UL_intro_AF_13$mutation_id <- paste0(UL_intro_AF_13$chr,"_",UL_intro_AF_13$pos)
UL_intro_AF_13 <- UL_intro_AF_13[grep("chr",UL_intro_AF_13$chr),]

UL_intro_AF_18 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/UL_2018.frq.txt")))
colnames(UL_intro_AF_18) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
UL_intro_AF_18$mutation_id <- paste0(UL_intro_AF_18$chr,"_",UL_intro_AF_18$pos)
UL_intro_AF_18 <- UL_intro_AF_18[grep("chr",UL_intro_AF_18$chr),]

UL_intro_tmp <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, p1 = GH_AF$REF_FREQ,
                           pmid =  UL_intro_AF_13$REF_FREQ, pt = UL_intro_AF_18$REF_FREQ)

head(UL_intro_tmp)

#get allele frequency change
UL_intro_tmp$sel_coef_1 <- (UL_intro_tmp$pmid - UL_intro_tmp$p1)
UL_intro_tmp$sel_coef_2 <- (UL_intro_tmp$pt - UL_intro_tmp$pmid)
UL_intro_tmp$sel_coef_total <- (UL_intro_tmp$pt - UL_intro_tmp$p1)

head(sel1)
#get only windows in that pop
UL_sel1 <- subset(sel1, in_UL == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
UL_intro_tmp$in_sel1 <- mapply(check_in_sel, UL_intro_tmp$chr, UL_intro_tmp$pos, MoreArgs = list(sel_df = UL_sel1))

#then just subset the alleles that are under selection
#update to include cut-offs from slims
UL_intro_tmp$in_sel1_allele <- ifelse(UL_intro_tmp$sel_coef_1 > 0.22 | UL_intro_tmp$sel_coef_1 < -0.21, TRUE, FALSE)
UL_intro_tmp$both_in_sel1 <- UL_intro_tmp$in_sel1 & UL_intro_tmp$in_sel1_allele

UL_sel_1_test <- subset(UL_intro_tmp, both_in_sel1 == "TRUE")

nrow(CA_sel_1_test)
nrow(UL_sel_1_test)

UL_sel1_plot <- ggplot(UL_sel_1_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  #  scale_fill_viridis_d(option = "C", direction = -1, na.value = "white") +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in the 1st Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)


UL_sel1_plot

#Now do sel2
#get only windows in that pop
UL_sel2 <- subset(sel2, in_UL == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
UL_intro_tmp$in_sel2 <- mapply(check_in_sel, UL_intro_tmp$chr, UL_intro_tmp$pos, MoreArgs = list(sel_df = UL_sel2))
head(UL_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
UL_intro_tmp$in_sel2_allele <- ifelse(UL_intro_tmp$sel_coef_2 > 0.20 | UL_intro_tmp$sel_coef_2 < -0.18, TRUE, FALSE)
UL_intro_tmp$both_in_sel2 <- UL_intro_tmp$in_sel2 & UL_intro_tmp$in_sel2_allele

UL_sel_2_test <- subset(UL_intro_tmp, both_in_sel2 == "TRUE")

nrow(CA_sel_2_test)
nrow(UL_sel_2_test)

UL_sel2_plot <- ggplot(UL_sel_2_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  #  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in 2nd Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)



UL_sel2_plot

#now total plot
#get only windows in that pop
UL_seltot <- subset(seltot, in_UL == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
UL_intro_tmp$in_seltot <- mapply(check_in_sel, UL_intro_tmp$chr, UL_intro_tmp$pos, MoreArgs = list(sel_df = UL_seltot))
head(UL_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
UL_intro_tmp$in_seltot_allele <- ifelse(UL_intro_tmp$sel_coef_total > 0.22 | UL_intro_tmp$sel_coef_total < -0.22, TRUE, FALSE)
UL_intro_tmp$both_in_seltot <- UL_intro_tmp$in_seltot & UL_intro_tmp$in_seltot_allele

UL_sel_tot_test <- subset(UL_intro_tmp, both_in_seltot == "TRUE")

nrow(CA_sel_tot_test)
nrow(UL_sel_tot_test)

UL_seltot_plot <- ggplot(UL_sel_tot_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  #  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in Total Time",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)



UL_seltot_plot

UL_combined_plot <- plot_grid(UL_sel1_plot, UL_sel2_plot, UL_seltot_plot, ncol=3, labels = c('(A)', '(B)', '(C)'))
UL_combined_plot

UL_final_plot <- ggdraw() +
  draw_label("UL", fontface = 'bold', size = 16, x = 0.5, y = 0.97, hjust = 0.5) +
  draw_plot(UL_combined_plot, y = 0, height = 0.95)

png("Documents/FIBR_temporal/allele_trajectory/UL_AF_density_plot.png", width = 10, height = 4, units = "in", res = 300)
UL_final_plot
dev.off()

#######now do it for Lower Laj
LL_intro_AF_13 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/LL_2013.frq.txt")))
colnames(LL_intro_AF_13) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
LL_intro_AF_13$mutation_id <- paste0(LL_intro_AF_13$chr,"_",LL_intro_AF_13$pos)
LL_intro_AF_13 <- LL_intro_AF_13[grep("chr",LL_intro_AF_13$chr),]

LL_intro_AF_18 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/LL_2018.frq.txt")))
colnames(LL_intro_AF_18) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
LL_intro_AF_18$mutation_id <- paste0(LL_intro_AF_18$chr,"_",LL_intro_AF_18$pos)
LL_intro_AF_18 <- LL_intro_AF_18[grep("chr",LL_intro_AF_18$chr),]

#head(GH_AF)
LL_intro_tmp <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, p1 = GH_AF$REF_FREQ,
                           pmid =  LL_intro_AF_13$REF_FREQ, pt = LL_intro_AF_18$REF_FREQ)


#get allele frequency change
LL_intro_tmp$sel_coef_1 <- (LL_intro_tmp$pmid - LL_intro_tmp$p1)
LL_intro_tmp$sel_coef_2 <- (LL_intro_tmp$pt - LL_intro_tmp$pmid)
LL_intro_tmp$sel_coef_total <- (LL_intro_tmp$pt - LL_intro_tmp$p1)

head(sel1)
#get only windows in that pop
LL_sel1 <- subset(sel1, in_LL == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
LL_intro_tmp$in_sel1 <- mapply(check_in_sel, LL_intro_tmp$chr, LL_intro_tmp$pos, MoreArgs = list(sel_df = LL_sel1))

#then just subset the alleles that are under selection
#update to include cut-offs from slims
LL_intro_tmp$in_sel1_allele <- ifelse(LL_intro_tmp$sel_coef_1 > 0.22 | LL_intro_tmp$sel_coef_1 < -0.22, TRUE, FALSE)
LL_intro_tmp$both_in_sel1 <- LL_intro_tmp$in_sel1 & LL_intro_tmp$in_sel1_allele

LL_sel_1_test <- subset(LL_intro_tmp, both_in_sel1 == "TRUE")

nrow(CA_sel_1_test)
nrow(LL_sel_1_test)

LL_sel1_plot <- ggplot(LL_sel_1_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  #  scale_fill_viridis_d(option = "C", direction = -1, na.value = "white") +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in 1st Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)


LL_sel1_plot

#Now do sel2
#get only windows in that pop
LL_sel2 <- subset(sel2, in_LL == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
LL_intro_tmp$in_sel2 <- mapply(check_in_sel, LL_intro_tmp$chr, LL_intro_tmp$pos, MoreArgs = list(sel_df = LL_sel2))
head(LL_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
LL_intro_tmp$in_sel2_allele <- ifelse(LL_intro_tmp$sel_coef_2 > 0.20 | LL_intro_tmp$sel_coef_2 < -0.18, TRUE, FALSE)
LL_intro_tmp$both_in_sel2 <- LL_intro_tmp$in_sel2 & LL_intro_tmp$in_sel2_allele

LL_sel_2_test <- subset(LL_intro_tmp, both_in_sel2 == "TRUE")

nrow(CA_sel_2_test)
nrow(LL_sel_2_test)

LL_sel2_plot <- ggplot(LL_sel_2_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  #  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in 2nd Time Interval",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)



LL_sel2_plot

#now total plot
#get only windows in that pop
LL_seltot <- subset(seltot, in_LL == "TRUE")

#code snps for whether they fall into a candidate loci and is represented by that pop
LL_intro_tmp$in_seltot <- mapply(check_in_sel, LL_intro_tmp$chr, LL_intro_tmp$pos, MoreArgs = list(sel_df = LL_seltot))
head(LL_intro_tmp)
#then just subset the alleles that are under selection
#update to include cut-offs from slims
LL_intro_tmp$in_seltot_allele <- ifelse(LL_intro_tmp$sel_coef_total > 0.23 | LL_intro_tmp$sel_coef_total < -0.22, TRUE, FALSE)
LL_intro_tmp$both_in_seltot <- LL_intro_tmp$in_seltot & LL_intro_tmp$in_seltot_allele

LL_sel_tot_test <- subset(LL_intro_tmp, both_in_seltot == "TRUE")

nrow(CA_sel_tot_test)
nrow(LL_sel_tot_test)

LL_seltot_plot <- ggplot(LL_sel_tot_test, aes(x = sel_coef_1, y = sel_coef_2)) +
  #  geom_point(color = "gray90", size = 0.3) +
  geom_density_2d_filled(contour_var = "density", bins = 15, alpha = 0.9) +
  scale_fill_viridis_d(option = "C", na.value = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_minimal()+
  labs(
    title = "Candidates in Total Time",
    x = "Delta AF in Time Period 1",
    y = "Delta AF in Time Period 2"
  ) +
  theme(legend.position = "none")+
  xlim(-0.8, 0.8) +
  ylim(-0.8, 0.8)



LL_seltot_plot

LL_combined_plot <- plot_grid(LL_sel1_plot, LL_sel2_plot, LL_seltot_plot, ncol=3, labels=c('(A)','(B)','(C)'))
LL_combined_plot

LL_final_plot <- ggdraw() +
  draw_label("LL", fontface = 'bold', size = 16, x = 0.5, y = 0.97, hjust = 0.5) +
  draw_plot(LL_combined_plot, y = 0, height = 0.95)

png("Documents/FIBR_temporal/allele_trajectory/LL_AF_density_plot.png", width = 10, height = 4, units = "in", res = 300)
LL_final_plot
dev.off()

LL_combined_plot
