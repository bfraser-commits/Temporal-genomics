library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(cowplot)
library(parallel)



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
CA_intro_tmp$p1_new <- ifelse(CA_intro_tmp$p1 > 0.5, 1 - CA_intro_tmp$p1, CA_intro_tmp$p1)
CA_intro_tmp$pmid_new <- ifelse(CA_intro_tmp$p1 > 0.5, 1 - CA_intro_tmp$pmid, CA_intro_tmp$pmid)
CA_intro_tmp$pt_new <- ifelse(CA_intro_tmp$p1 > 0.5, 1 - CA_intro_tmp$pt, CA_intro_tmp$pt)

#get allele frequency change
#CA_intro_tmp$sel_coef_1 <- (CA_intro_tmp$pmid - CA_intro_tmp$p1)
#CA_intro_tmp$sel_coef_2 <- (CA_intro_tmp$pt - CA_intro_tmp$pmid)
#CA_intro_tmp$sel_coef_total <- (CA_intro_tmp$pt - CA_intro_tmp$p1)

CA_intro_tmp$sel_coef_1 <- (CA_intro_tmp$pmid_new - CA_intro_tmp$p1_new)
CA_intro_tmp$sel_coef_2 <- (CA_intro_tmp$pt_new - CA_intro_tmp$pmid_new)
CA_intro_tmp$sel_coef_total <- (CA_intro_tmp$pt_new - CA_intro_tmp$p1_new)


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
#TA_intro_tmp$sel_coef_1 <- (TA_intro_tmp$pmid - TA_intro_tmp$p1)
#TA_intro_tmp$sel_coef_2 <- (TA_intro_tmp$pt - TA_intro_tmp$pmid)
#TA_intro_tmp$sel_coef_total <- (TA_intro_tmp$pt - TA_intro_tmp$p1)

#switch to minor allele
TA_intro_tmp$p1_new <- ifelse(TA_intro_tmp$p1 > 0.5, 1 - TA_intro_tmp$p1, TA_intro_tmp$p1)
TA_intro_tmp$pmid_new <- ifelse(TA_intro_tmp$p1 > 0.5, 1 - TA_intro_tmp$pmid, TA_intro_tmp$pmid)
TA_intro_tmp$pt_new <- ifelse(TA_intro_tmp$p1 > 0.5, 1 - TA_intro_tmp$pt, TA_intro_tmp$pt)

TA_intro_tmp$sel_coef_1 <- (TA_intro_tmp$pmid_new - TA_intro_tmp$p1_new)
TA_intro_tmp$sel_coef_2 <- (TA_intro_tmp$pt_new - TA_intro_tmp$pmid_new)
TA_intro_tmp$sel_coef_total <- (TA_intro_tmp$pt_new - TA_intro_tmp$p1_new)




#UL
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
#UL_intro_tmp$sel_coef_1 <- (UL_intro_tmp$pmid - UL_intro_tmp$p1)
#UL_intro_tmp$sel_coef_2 <- (UL_intro_tmp$pt - UL_intro_tmp$pmid)
#UL_intro_tmp$sel_coef_total <- (UL_intro_tmp$pt - UL_intro_tmp$p1)

#switch to minor allele
UL_intro_tmp$p1_new <- ifelse(UL_intro_tmp$p1 > 0.5, 1 - UL_intro_tmp$p1, UL_intro_tmp$p1)
UL_intro_tmp$pmid_new <- ifelse(UL_intro_tmp$p1 > 0.5, 1 - UL_intro_tmp$pmid, UL_intro_tmp$pmid)
UL_intro_tmp$pt_new <- ifelse(UL_intro_tmp$p1 > 0.5, 1 - UL_intro_tmp$pt, UL_intro_tmp$pt)

UL_intro_tmp$sel_coef_1 <- (UL_intro_tmp$pmid_new - UL_intro_tmp$p1_new)
UL_intro_tmp$sel_coef_2 <- (UL_intro_tmp$pt_new - UL_intro_tmp$pmid_new)
UL_intro_tmp$sel_coef_total <- (UL_intro_tmp$pt_new - UL_intro_tmp$p1_new)

#LL
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
#LL_intro_tmp$sel_coef_1 <- (LL_intro_tmp$pmid - LL_intro_tmp$p1)
#LL_intro_tmp$sel_coef_2 <- (LL_intro_tmp$pt - LL_intro_tmp$pmid)
#LL_intro_tmp$sel_coef_total <- (LL_intro_tmp$pt - LL_intro_tmp$p1)

#switch to minor allele
LL_intro_tmp$p1_new <- ifelse(LL_intro_tmp$p1 > 0.5, 1 - LL_intro_tmp$p1, LL_intro_tmp$p1)
LL_intro_tmp$pmid_new <- ifelse(LL_intro_tmp$p1 > 0.5, 1 - LL_intro_tmp$pmid, LL_intro_tmp$pmid)
LL_intro_tmp$pt_new <- ifelse(LL_intro_tmp$p1 > 0.5, 1 - LL_intro_tmp$pt, LL_intro_tmp$pt)

LL_intro_tmp$sel_coef_1 <- (LL_intro_tmp$pmid_new - LL_intro_tmp$p1_new)
LL_intro_tmp$sel_coef_2 <- (LL_intro_tmp$pt_new - LL_intro_tmp$pmid_new)
LL_intro_tmp$sel_coef_total <- (LL_intro_tmp$pt_new - LL_intro_tmp$p1_new)



head(TA_intro_tmp)
####merge and get sum change

Total_dataset_sel1 <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, CA_sel1 = CA_intro_tmp$sel_coef_1,
                            TA_sel1 = TA_intro_tmp$sel_coef_1, LL_sel1 = LL_intro_tmp$sel_coef_1, UL_sel1 = UL_intro_tmp$sel_coef_1)
summary(Total_dataset_sel1$sum_sel1)

#Total_dataset_sel1$sum_sel1 <- (abs(Total_dataset_sel1$CA_sel1) + abs(Total_dataset_sel1$TA_sel1) + abs(Total_dataset_sel1$LL_sel1) + abs(Total_dataset_sel1$UL_sel1))
Total_dataset_sel1$sum_sel1 <- abs(Total_dataset_sel1$CA_sel1 + Total_dataset_sel1$TA_sel1 + Total_dataset_sel1$LL_sel1 + Total_dataset_sel1$UL_sel1)

Total_dataset_sel2 <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, CA_sel2 = CA_intro_tmp$sel_coef_2,
                                 TA_sel2 = TA_intro_tmp$sel_coef_2, LL_sel2 = LL_intro_tmp$sel_coef_2, UL_sel2 = UL_intro_tmp$sel_coef_2)

#Total_dataset_sel2$sum_sel2 <- (abs(Total_dataset_sel2$CA_sel2) + abs(Total_dataset_sel2$TA_sel2) + abs(Total_dataset_sel2$LL_sel2) + abs(Total_dataset_sel2$UL_sel2))
Total_dataset_sel2$sum_sel2 <- abs(Total_dataset_sel2$CA_sel2 + Total_dataset_sel2$TA_sel2 + Total_dataset_sel2$LL_sel2 + Total_dataset_sel2$UL_sel2)
summary(Total_dataset_sel2$sum_sel2)

Total_dataset_seltot <- data.frame(chr = GH_AF$chr, pos = GH_AF$pos, mutation_id = GH_AF$mutation_id, CA_seltot = CA_intro_tmp$sel_coef_total,
                                 TA_seltot = TA_intro_tmp$sel_coef_total, LL_seltot = LL_intro_tmp$sel_coef_total, UL_seltot = UL_intro_tmp$sel_coef_total)

#Total_dataset_seltot$sum_seltot <- (abs(Total_dataset_seltot$CA_seltot) + abs(Total_dataset_seltot$TA_seltot) + abs(Total_dataset_seltot$LL_seltot) + abs(Total_dataset_seltot$UL_seltot))
Total_dataset_seltot$sum_seltot <- abs(Total_dataset_seltot$CA_seltot + Total_dataset_seltot$TA_seltot + Total_dataset_seltot$LL_seltot + Total_dataset_seltot$UL_seltot)

summary(Total_dataset_sel1$sum_sel1)
summary(Total_dataset_sel2$sum_sel2)
summary(Total_dataset_seltot$sum_seltot)

##get em-pvalues
slims_sel1 <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/summary_sel_coef1.txt", header=T)
simulated_max <- slims_sel1$quantile_95

Total_dataset_sel1$empirical_p <- sapply(Total_dataset_sel1$sum_sel1,  function(x) {mean(simulated_max >= x, na.rm = TRUE) })
out1 <- subset(Total_dataset_sel1, empirical_p < 0.05)
nrow(out1)
nrow(out1)/nrow(Total_dataset_sel1)


sel1_emp_p <- ggplot(Total_dataset_sel1, aes(x = empirical_p)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  labs(
    title = expression(1^st ~ "Time Interval"),
    x = "Empirical P-value",
    y = "Density"
  )+
  theme_bw()

sel1_emp_p

#selcoef2
slims_sel2 <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/summary_sel_coef2.txt", header=T)
simulated_max2 <- slims_sel2$quantile_95

Total_dataset_sel2$empirical_p <- sapply(Total_dataset_sel2$sum_sel2,  function(x) {mean(simulated_max2 >= x, na.rm = TRUE) })

out2 <- subset(Total_dataset_sel2, empirical_p < 0.05)
nrow(out2)
nrow(out2)/nrow(Total_dataset_sel2)

sel2_emp_p <- ggplot(Total_dataset_sel2, aes(x = empirical_p)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  labs(
    title = expression(2^nd ~ "Time Interval"),
    x = "Empirical P-value",
    y = "Density"
  )+
  theme_bw()

sel2_emp_p

#seltot
slims_seltot <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/summary_sel_tot.txt", header=T)
simulated_maxtot <- slims_seltot$quantile_95

Total_dataset_seltot$empirical_p <- sapply(Total_dataset_seltot$sum_seltot,  function(x) {mean(simulated_maxtot >= x, na.rm = TRUE) })

outtot <- subset(Total_dataset_seltot, empirical_p < 0.05)
nrow(outtot)
nrow(outtot)/nrow(Total_dataset_seltot)

seltot_emp_p <-  ggplot(Total_dataset_seltot, aes(x = empirical_p)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  labs(
    title = "Total Time Interval",
    x = "Empirical P-value",
    y = "Density"
  )+
  theme_bw()

png(filename = "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/emp_pvalue.png",  width = 8, height = 12, units = "in",   res = 300)
plot_grid(sel1_emp_p, sel2_emp_p, seltot_emp_p, ncol=1)
dev.off()

#output final figs
sel1_p <-  ggplot(Total_dataset_sel1, aes(x = sum_sel1)) +
  geom_density(fill = "steelblue", alpha = 0.5, adjust = 1.5) +
  labs(
    title = expression(1^st ~ "Time Interval"),
    x = "summed allele frequency change",
    y = "Density"
  )+
  theme_bw()

sel1_p

sel2_p <-  ggplot(Total_dataset_sel2, aes(x = sum_sel2)) +
  geom_density(fill = "steelblue", alpha = 0.5, adjust = 1.5) +
  labs(
    title = expression(2^nd ~ "Time Interval"),
    x = "summed allele frequency change",
    y = "Density"
  )+
  theme_bw()

sel2_p

seltot_p <-  ggplot(Total_dataset_seltot, aes(x = sum_seltot)) +
  geom_density(fill = "steelblue", alpha = 0.5, adjust = 1.5) +
  labs(
    title = "Total Time Interval",
    x = "summed allele frequency change",
    y = "Density"
  )+
  theme_bw()

seltot_p

png(filename = "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/summed_change.png",  width = 8, height = 12, units = "in",   res = 300)
plot_grid(sel1_p, sel2_p, seltot_p, ncol=1)
dev.off()

###windowise p-values
## average the per SNP values into windows
head(Total_dataset_sel1)

chrs<-unique(Total_dataset_sel1$chr)
wind_size<-50000

winds<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-Total_dataset_sel1[Total_dataset_sel1$chr == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_p<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$pos <= winds2[y] & tmp$pos >= winds1[y],]
    #out<-data.frame(mean(tmp2$empirical_p,na.rm = T))
    out <- data.frame(sum(tmp2$empirical_p < 0.05, na.rm = TRUE))
    # Tidy
    out$total<-length(tmp2$pos)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    #out$comp<-'UL18'
    colnames(out)<-c('pvalue_count','total_size','chrom','window','BP1','BP2')
    
    return(out)
  })))
  return(sum_p)
},mc.cores=3)))


#do it for sel2
chrs<-unique(Total_dataset_sel2$chr)
wind_size<-50000

winds_sel2<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-Total_dataset_sel2[Total_dataset_sel2$chr == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_p<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$pos <= winds2[y] & tmp$pos >= winds1[y],]
    #out<-data.frame(mean(tmp2$empirical_p,na.rm = T))
    out <- data.frame(sum(tmp2$empirical_p < 0.05, na.rm = TRUE))
    # Tidy
    out$total<-length(tmp2$pos)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    #out$comp<-'UL18'
    colnames(out)<-c('pvalue_count','total_size','chrom','window','BP1','BP2')
    
    return(out)
  })))
  return(sum_p)
},mc.cores=3)))

###seltot
chrs<-unique(Total_dataset_seltot$chr)
wind_size<-50000

winds_seltot<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-Total_dataset_seltot[Total_dataset_seltot$chr == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_p<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$pos <= winds2[y] & tmp$pos >= winds1[y],]
    #out<-data.frame(mean(tmp2$empirical_p,na.rm = T))
    out <- data.frame(sum(tmp2$empirical_p < 0.05, na.rm = TRUE))
    # Tidy
    out$total<-length(tmp2$pos)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    #out$comp<-'UL18'
    colnames(out)<-c('pvalue_count','total_size','chrom','window','BP1','BP2')
    
    return(out)
  })))
  return(sum_p)
},mc.cores=3)))



#calculate outlier windows
winds$fold_enrichment <-  (winds$pvalue_count / winds$total_size) / 0.05
winds_test <- subset(winds, total_size < 1000 & total_size > 50)
winds_test_out <- subset(winds_test, fold_enrichment > 5)
nrow(winds_test_out)
table(winds_test_out$chrom)
#subset(winds_test_out, chrom == "chr15")

winds_sel2$fold_enrichment <- (winds_sel2$pvalue_count / winds_sel2$total_size) / 0.05
winds_test_sel2 <- subset(winds_sel2, total_size < 1000 & total_size > 50)
winds_test_sel2_out <- subset(winds_test_sel2, fold_enrichment > 5)
nrow(winds_test_sel2_out)
table(winds_test_sel2_out$chrom)

winds_seltot$fold_enrichment <- (winds_seltot$pvalue_count / winds_seltot$total_size) / 0.05
winds_test_seltot <- subset(winds_seltot, total_size < 1000 & total_size > 50)
winds_test_seltot_out <- subset(winds_test_seltot, fold_enrichment > 5)
nrow(winds_test_seltot_out)
table(winds_test_seltot_out$chrom)

#check for overlap among time interval outs
winds_test_out$ID <- paste(winds_test_out$chrom, winds_test_out$BP1, sep = "_")
head(winds_test_out)
winds_test_sel2_out$ID <- paste(winds_test_sel2_out$chrom, winds_test_sel2_out$BP1, sep = "_")
head(winds_test_sel2_out)
winds_test_seltot_out$ID <- paste(winds_test_seltot_out$chrom, winds_test_seltot_out$BP1, sep = "_")
head(winds_test_seltot_out)
#merge
merge_1_2 <- merge(winds_test_out, winds_test_sel2_out, by = "ID")
nrow(merge_1_2)
merge_1_tot <- merge(winds_test_out, winds_test_seltot_out, by = "ID")
nrow(merge_1_tot)
merge_2_tot <- merge(winds_test_sel2_out, winds_test_seltot_out, by = "ID")
nrow(merge_2_tot)


table(winds_test_out$chrom)
table(winds_test_sel2_out$chrom)
table(winds_test_seltot_out$chrom)

#write outputs
head(Total_dataset_sel1)

write.table(Total_dataset_sel1, file = "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/Total_dataset_sel1.txt", sep = "\t", quote=F, row.names=F)
write.table(Total_dataset_sel2, file = "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/Total_dataset_sel2.txt", sep = "\t", quote=F, row.names=F)
write.table(Total_dataset_seltot, file = "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/Total_dataset_seltot.txt", sep = "\t", quote=F, row.names=F)

#write windows
write.table(winds, file= "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/winds_sel1.txt", sep = "\t", quote=F, row.names=F)
write.table(winds_sel2, file= "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/winds_sel2.txt", sep = "\t", quote=F, row.names=F)
write.table(winds_seltot, file= "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/winds_seltot.txt", sep = "\t", quote=F, row.names=F)

#
winds <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/winds_sel1.txt", header=T)
winds_sel2 <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/winds_sel2.txt", header=T)
winds_seltot <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/winds_seltot.txt", header=T)
head(winds)
###full genome plots###
#sel1
# Order chromosomes numerically (handles names like chr1 or 1)
chr_names <- unique(as.character(winds_test$chrom))
chr_num <- as.numeric(sub("^chr", "", chr_names))
chr_order <- chr_names[order(chr_num)]

# Make chromosome an ordered factor
winds_test$chrom <- factor(as.character(winds_test$chrom),  levels = chr_order)

# Chromosome lengths in correct order
chr_len <- tapply(winds_test$BP1, winds_test$chrom, max)
chr_len <- chr_len[chr_order]

# Cumulative offsets
offset <- c(0, cumsum(chr_len[-length(chr_len)]))
names(offset) <- names(chr_len)

# Add cumulative positions
winds_test$BPcum <- winds_test$BP1 +  offset[as.character(winds_test$chrom)]

# Chromosome centres for axis labels
axisdf <- aggregate(BPcum ~ chrom,  winds_test,  function(x) (min(x) + max(x)) / 2)

axisdf <- axisdf[match(chr_order, as.character(axisdf$chrom)), ]

# Alternating grey colours
chrom_cols <- setNames( rep(c("grey70", "grey30"), length.out = length(chr_order)),  chr_order)

# Plot
wg_plot <- ggplot(winds_test, aes(BPcum, fold_enrichment)) +
  geom_point(aes(colour = ifelse(fold_enrichment > 5, "red",as.character(chrom))),size = 0.6) +
  scale_colour_manual(values = c(chrom_cols, "red" = "red"),guide = "none") +
  scale_x_continuous(breaks = axisdf$BPcum,labels = axisdf$chrom) +
  labs(title = expression("(A) " * 1^st ~ "Period"),
    x = "Chromosome",
    y = "Fold enrichment"
  ) +
  coord_cartesian(ylim = c(0, 11.01)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

wg_plot


#sel2
# Order chromosomes numerically (handles names like chr1 or 1)
chr_names <- unique(as.character(winds_test_sel2$chrom))
chr_num <- as.numeric(sub("^chr", "", chr_names))
chr_order <- chr_names[order(chr_num)]

# Make chromosome an ordered factor
winds_test_sel2$chrom <- factor(
  as.character(winds_test_sel2$chrom),
  levels = chr_order
)

# Chromosome lengths in correct order
chr_len <- tapply(winds_test_sel2$BP1, winds_test_sel2$chrom, max)
chr_len <- chr_len[chr_order]

# Cumulative offsets
offset <- c(0, cumsum(chr_len[-length(chr_len)]))
names(offset) <- names(chr_len)

# Add cumulative positions
winds_test_sel2$BPcum <- winds_test_sel2$BP1 +
  offset[as.character(winds_test_sel2$chrom)]

# Chromosome centres for axis labels
axisdf <- aggregate(
  BPcum ~ chrom,
  winds_test_sel2,
  function(x) (min(x) + max(x)) / 2
)

axisdf <- axisdf[match(chr_order, as.character(axisdf$chrom)), ]

# Alternating grey colours
chrom_cols <- setNames(
  rep(c("grey70", "grey30"), length.out = length(chr_order)),
  chr_order
)

# Plot
wg_plot_sel2 <- ggplot(winds_test_sel2, aes(BPcum, fold_enrichment)) +
  geom_point(
    aes(colour = ifelse(
      fold_enrichment > 5,
      "red",
      as.character(chrom)
    )),
    size = 0.6
  ) +
  scale_colour_manual(
    values = c(chrom_cols, "red" = "red"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = axisdf$BPcum,
    labels = axisdf$chrom
  ) +
  labs(title = expression("(B) " *2^nd ~ "Period"),
    x = "Chromosome",
    y = "Fold enrichment"
  ) +
  coord_cartesian(ylim = c(0, 11.01)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

wg_plot_sel2

#seltot
# Order chromosomes numerically (handles names like chr1 or 1)
chr_names <- unique(as.character(winds_test_seltot$chrom))
chr_num <- as.numeric(sub("^chr", "", chr_names))
chr_order <- chr_names[order(chr_num)]

# Make chromosome an ordered factor
winds_test_seltot$chrom <- factor(
  as.character(winds_test_seltot$chrom),
  levels = chr_order
)

# Chromosome lengths in correct order
chr_len <- tapply(winds_test_seltot$BP1, winds_test_seltot$chrom, max)
chr_len <- chr_len[chr_order]

# Cumulative offsets
offset <- c(0, cumsum(chr_len[-length(chr_len)]))
names(offset) <- names(chr_len)

# Add cumulative positions
winds_test_seltot$BPcum <- winds_test_seltot$BP1 +
  offset[as.character(winds_test_seltot$chrom)]

# Chromosome centres for axis labels
axisdf <- aggregate(
  BPcum ~ chrom,
  winds_test_seltot,
  function(x) (min(x) + max(x)) / 2
)

axisdf <- axisdf[match(chr_order, as.character(axisdf$chrom)), ]

# Alternating grey colours
chrom_cols <- setNames(
  rep(c("grey70", "grey30"), length.out = length(chr_order)),
  chr_order
)

# Plot
wg_plot_seltot <- ggplot(winds_test_seltot, aes(BPcum, fold_enrichment)) +
  geom_point(
    aes(colour = ifelse(
      fold_enrichment > 5,
      "red",
      as.character(chrom)
    )),
    size = 0.6
  ) +
  scale_colour_manual(
    values = c(chrom_cols, "red" = "red"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = axisdf$BPcum,
    labels = axisdf$chrom
  ) +
  labs(title = "(C) Total Time",
    x = "Chromosome",
    y = "Fold enrichment"
  ) +
  coord_cartesian(ylim = c(0, 11.01)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

wg_plot_seltot
max(winds_test$fold_enrichment)
max(winds_test_sel2$fold_enrichment)

png("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/WG_plots.png", height=6, width=8, units = "in",   res = 300)
plot_grid(wg_plot, wg_plot_sel2, wg_plot_seltot, ncol=1)
dev.off()

###chrom plots
chr_sel1 <- subset(winds_test, chrom == "chr20")
chr_sel2 <- subset(winds_test_sel2, chrom == "chr20")
chr_seltot <- subset(winds_test_seltot, chrom == "chr20")


chrom_sel1_plot <- ggplot(chr_sel1, aes(x=(BP1+500), y=fold_enrichment))+
  geom_point(aes(color = fold_enrichment > 5 )) +
  scale_color_manual(
    values = c("FALSE" = "grey", "TRUE" = "red")
  ) +
  labs(title = expression("(D) chr20 - " * 1^st ~ "Period"),
       x = "Chromosome",
       y = "Fold enrichment"
  ) +
  coord_cartesian(ylim = c(0, 11.01)) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

chrom_sel1_plot

chrom_sel2_plot <- ggplot(chr_sel2, aes(x=(BP1+500), y=fold_enrichment))+
  geom_point(aes(color = fold_enrichment > 5 )) +
  scale_color_manual(
    values = c("FALSE" = "grey", "TRUE" = "red")
  ) +
  labs(title = expression("(E) chr20 - " * 2^nd ~ "Period"),
       x = "Chromosome",
       y = "Fold enrichment"
  ) +
  coord_cartesian(ylim = c(0, 11.01)) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )


chrom_sel2_plot

chrom_seltot_plot <- ggplot(chr_seltot, aes(x=(BP1+500), y=fold_enrichment))+
  geom_point(aes(color = fold_enrichment > 5 )) +
  scale_color_manual(
    values = c("FALSE" = "grey", "TRUE" = "red")
  ) +
  labs(title = "(F) chr20 - Total Time",
       x = "Chromosome",
       y = "Fold enrichment"
  ) +
  coord_cartesian(ylim = c(0, 11.01)) +
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )


chrom_seltot_plot
png("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/chr20_plots.png", height=6, width=4, units = "in",   res = 300)
plot_grid(chrom_sel1_plot, chrom_sel2_plot, chrom_seltot_plot, ncol=1)
dev.off()
