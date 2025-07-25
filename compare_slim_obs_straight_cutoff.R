slim <- data.frame(fread("Documents/FIBR_temporal/slims_test/all_pop_sims_min_13_18.txt", header=T))
I_slim <- subset(slim, pop == "UL")
I_slim$sel_coef_1 <- (I_slim$pmid - I_slim$p1)
I_slim$sel_coef_2 <- (I_slim$pt - I_slim$pmid)
I_slim$sel_coef_total <- (I_slim$pt - I_slim$p1)

I_slim <- I_slim[!(I_slim$p1 == 1 & I_slim$pt == 1), ]
#I_slim <- I_slim[I_slim$p1 >= 0.05 & I_slim$p1 <= 0.95, ]

quantile(abs(I_slim$sel_coef_total), 0.99, na.rm = TRUE)
quantile(abs(I_slim$sel_coef_1), 0.99, na.rm = TRUE)
quantile(abs(I_slim$sel_coef_2), 0.99, na.rm = TRUE)


# Per pop S
GH_AF <- data.frame(fread("Documents/FIBR_temporal/slims_test/GH.frq.txt"))
#GH_AF <- data.frame(fread("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/AF_freqs/GH.frq.txt"))
colnames(GH_AF) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
GH_AF$mutation_id <- paste0(GH_AF$chr,"_",GH_AF$pos)
GH_AF <- GH_AF[grep("chr",GH_AF$chr),]

# Fetch the minor af as lowest row min of both freqs...
#GH_maf <- rowMins(as.matrix(GH_AF[,c("REF_FREQ","ALT_FREQ")]),value=T)
#GH_maf_index <- rowMins(as.matrix(GH_AF[,c("REF_FREQ","ALT_FREQ")]),value=F)

intro_AF_13 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/UL_2013.frq.txt")))
colnames(intro_AF_13) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
intro_AF_13$mutation_id <- paste0(intro_AF_13$chr,"_",intro_AF_13$pos)
intro_AF_13 <- intro_AF_13[grep("chr",intro_AF_13$chr),]

intro_AF_18 <- data.frame(fread(paste0("Documents/FIBR_temporal/slims_test/UL_2018.frq.txt")))
colnames(intro_AF_18) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
intro_AF_18$mutation_id <- paste0(intro_AF_18$chr,"_",intro_AF_18$pos)
intro_AF_18 <- intro_AF_18[grep("chr",intro_AF_18$chr),]

intro_tmp <- data.frame(mutation_id = GH_AF$mutation_id, p1 = GH_AF$REF_FREQ,
                        pmid =  intro_AF_13$REF_FREQ, pt = intro_AF_18$REF_FREQ)


intro_tmp$sel_coef_1 <- (intro_tmp$pmid - intro_tmp$p1)
intro_tmp$sel_coef_2 <- (intro_tmp$pt - intro_tmp$p1)
intro_tmp$sel_coef_total <- (intro_tmp$pt - intro_tmp$p1)

summary(intro_tmp)
summary(I_slim)

#remove any sites that were fixed in GH

#intro_tmp <- subset(intro_tmp, p1 != 0)

#remove any sites that were fixed at the beginning and end
intro_tmp <- intro_tmp[!( (intro_tmp$p1 == 0 & intro_tmp$pt == 0) | 
                            (intro_tmp$p1 == 1 & intro_tmp$pt == 1) ), ]

#intro_tmp <- intro_tmp[intro_tmp$p1 >= 0.05 & intro_tmp$p1 <= 0.95, ]



head(intro_tmp)
head(I_slim)

#now do the starting AF to final AF in slims 
I_slim <- I_slim %>%
  mutate(start_bin = cut(p1,
                         breaks = seq(0, 1, by = 0.1),
                         include.lowest = TRUE,
                         right = FALSE,
                         labels = paste0(seq(0, 0.9, by = 0.1), "-", seq(0.1, 1, by = 0.1))))


summary(I_slim)

intro_tmp <- intro_tmp %>%
  mutate(start_bin = cut(p1,
                         breaks = seq(0, 1, by = 0.1),
                         include.lowest = TRUE,
                         right = FALSE,
                         labels = paste0(seq(0, 0.9, by = 0.1), "-", seq(0.1, 1, by = 0.1))))

intro_tmp$dataset <- "observed"
I_slim$dataset <- "simulated"
head(intro_tmp)
head(I_slim)
I_slim <- data.frame(I_slim$mutation_id, I_slim$p1, I_slim$pmid, I_slim$pt, I_slim$sel_coef_1, I_slim$sel_coef_2, I_slim$sel_coef_total, I_slim$start_bin, I_slim$dataset)
head(I_slim)
colnames(I_slim) <-c ("mutation_id", "p1", "pmid","pt", "sel_coef_1", "sel_coef_2", "sel_coef_total", "start_bin", "dataset")
final_df <- data.frame(rbind(intro_tmp, I_slim))
head(I_slim)

sup_fig <- ggplot(final_df, aes(x = dataset, y = sel_coef_1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~ start_bin, ncol = 3) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  labs(
    x = "Dataset",
    y = "Change in Allele Frequency (ΔAF)",
    fill = "Dataset",
    title = "Comparison of ΔAF Distributions by Starting Frequency Bin"
  ) +
  theme_minimal()


sup_fig

sup_fig <- ggplot(final_df, aes(x = dataset, y = sel_coef_1, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  #  facet_wrap(~ start_bin, ncol = 3) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  labs(
    x = "Dataset",
    y = "Change in Allele Frequency (ΔAF)",
    fill = "Dataset",
    title = "Comparison of ΔAF Distributions by Starting Frequency Bin"
  ) +
  theme_minimal()

sup_fig_2 <- ggplot(final_df, aes(x = dataset, y = sel_coef_2, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  #  facet_wrap(~ start_bin, ncol = 3) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  labs(
    x = "Dataset",
    y = "Change in Allele Frequency (ΔAF)",
    fill = "Dataset",
    title = "Comparison of ΔAF Distributions by Starting Frequency Bin"
  ) +
  theme_minimal()

sup_fig_2


sup_fig_tot <- ggplot(final_df, aes(x = dataset, y = sel_coef_total, fill = dataset)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
  #  facet_wrap(~ start_bin, ncol = 3) +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  labs(
    x = "Dataset",
    y = "Change in Allele Frequency (ΔAF)",
    fill = "Dataset",
    title = "Comparison of ΔAF Distributions by Starting Frequency Bin"
  ) +
  theme_minimal()

sup_fig_tot


pdf("Documents/FIBR_temporal/slims_test/LL_sup_figure.pdf", width = 10, heigh=8)
sup_fig
dev.off()

head(final_df)
low_af <- subset(final_df, start_bin == "0.4-0.5" )
summary(low_af)

ggplot(final_df, aes(x = sel_coef_1, fill = dataset)) +
  geom_density(alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("observed" = "steelblue", "simulated" = "orange")) +
  labs(
    x = "Change in Allele Frequency (ΔAF)",
    y = "Density",
    fill = "Dataset",
    title = "Comparison of ΔAF Distributions by Starting Frequency Bin"
  ) +
  theme_minimal()


# Can now take a cutoff from this distribution in absolute terms and highlight SNPs
#try with 99
sim_cutoff_sel1 <- quantile(abs(I_slim$sel_coef_1),0.99)
sim_cutoff_sel2 <- quantile(abs(I_slim$sel_coef_2),0.99)
sim_cutoff_seltot <- quantile(abs(I_slim$sel_coef_total),0.99)

sim_cutoff_sel1 
sim_cutoff_sel2 
sim_cutoff_seltot

#get outliers
split_vals <- do.call(rbind, strsplit(intro_tmp$mutation_id, "_"))
intro_tmp$chr <- split_vals[, 1]
intro_tmp$pos <- as.numeric(split_vals[, 2])

head(intro_tmp)
# Get outlier SNPs
#summary(intro_tmp)
outlier_snps_coef1 <- intro_tmp[abs(intro_tmp$sel_coef_1) > sim_cutoff_sel1,]
nrow(outlier_snps_coef1)
(nrow(outlier_snps_coef1)/nrow(intro_tmp))*100
outlier_snps_coef2 <- intro_tmp[abs(intro_tmp$sel_coef_2) > sim_cutoff_sel2,]
nrow(outlier_snps_coef2)
(nrow(outlier_snps_coef2)/nrow(intro_tmp))*100
outlier_snps_coeftot <- intro_tmp[abs(intro_tmp$sel_coef_total) > sim_cutoff_seltot,]
nrow(outlier_snps_coeftot)
(nrow(outlier_snps_coeftot)/nrow(intro_tmp))*100

# How are these organised among chromosomes
table(outlier_snps_coef1$chr)
table(outlier_snps_coef2$chr)
table(outlier_snps_coeftot$chr)

head(intro_tmp)

#windowise
chrs<-unique(intro_tmp$chr)
wind_size<-50000
#chrs

winds<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-intro_tmp[intro_tmp$chr == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$pos <= winds2[y] & tmp$pos >= winds1[y],]
    if (nrow(tmp2) > 0) {
      outlier_count_sel1 <- sum(abs(tmp2$sel_coef_1) > sim_cutoff_sel1, na.rm = TRUE) 
      outlier_count_sel2 <- sum(abs(tmp2$sel_coef_2) > sim_cutoff_sel2, na.rm = TRUE) 
      outlier_count_seltot <- sum(abs(tmp2$sel_coef_total) > sim_cutoff_seltot, na.rm = TRUE) 
      site_count <- length(tmp2$sel_coef_1)
      mean_sel1 <- mean(abs(tmp2$sel_coef_1)) 
      mean_sel2 <- mean(abs(tmp2$sel_coef_2)) 
      mean_seltot <- mean(abs(tmp2$sel_coef_total)) 
      out <- data.frame(outlier_count_sel1 = outlier_count_sel1, outlier_count_sel2 = outlier_count_sel2, outlier_count_seltot = outlier_count_seltot, site_count = site_count, mean_sel1 = mean_sel1, mean_sel2 = mean_sel2, mean_seltot = mean_seltot)
    } else {
      out <- data.frame(outlier_count_sel1 = NA, outlier_count_sel2 = NA, outlier_count_seltot = NA, site_count = NA, mean_sel1 = NA, mean_sel2 = NA, mean_seltot = NA)  # Handle empty windows
    }
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'UL'
    colnames(out)<-c('outlier_count_sel1','outlier_count_sel2', 'outlier_count_seltot', 'site_count', "mean_sel1", "mean_sel2", "mean_seltot", 'chr','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

head(winds)
summary(winds)
summary(winds$mean_sel2)
table(winds$chr)

write.table(winds, file = "Documents/FIBR_temporal/UL_99straight_cutoff_50kb.txt", quote=F, row.names=F)

###test
chr <- subset(winds, chr == "chr15")
head(chr)

chrom_plot <- ggplot(chr, aes(x=(BP1+500), y=outlier_count_sel1))+
  geom_point(size = 0.5)+
 geom_point(aes(x = BP1 + 500, y = outlier_count_sel2), color = "red", size = 0.5)+ 
  geom_point(aes(x = BP1 + 500, y = outlier_count_seltot), color = "blue", size = 0.5)+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()

chrom_plot

chrom_plot <- ggplot(chr, aes(x=(BP1+500), y=mean_sel1))+
  geom_point(size = 0.5)+
#  geom_point(aes(x = BP1 + 500, y = outlier_count_sel2), color = "red", size = 0.5)+ 
 # geom_point(aes(x = BP1 + 500, y = outlier_count_seltot), color = "blue", size = 0.5)+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()
