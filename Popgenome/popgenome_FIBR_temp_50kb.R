### Rscript for popGenome fst, pi and tajD ###
#install.packages("ff",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')
#install.packages("data.table",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')
#install.packages("dplyr",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')
#install.packages("PopGenome",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')

library("ff", lib.loc='~/bin/R/library')
library("data.table", lib.loc='~/bin/R/library')
library("PopGenome", lib.loc='~/bin/R/library')
#library("dplyr", lib.loc='~/bin/R/library')

##UPDATE WITH FAI FOR WHATEVER VCF YOURE USING
#Import fai for genome data
chr_length<-read.table("/gpfs/ts0/projects/Research_Project-T109423/STAR/STAR.chromosomes.release.fasta.fai",header=F)

####### Loop over chromosomes ###########
chr_function<-function(x){

#chr name as defined from fasta.fai
tid=as.character(chr_length[x,]$V1)


#chr length
topos = chr_length[x,]$V2


#Read in VCF
vcf_file <- readVCF('/gpfs/ts0/projects/Research_Project-T109423/ali/FIBR_transfer/ftp1.sequencing.exeter.ac.uk/V0085/11_trimmed/PopGenome/FIBR_2018_liftover_phased.vcf.gz', tid=tid, frompos = 1, topos = topos, numcols = 1000000, include.unknown=TRUE)

### assign pops ###
CA13 <- c('C12','C13','C15','C20','C21','C23','C4','C6','C8','C9')
CA18 <- c('CA10_18','CA11_18','CA12_18','CA13_18','CA14_18','CA15_18','CA16_18','CA17_18','CA18_18','CA19_18','CA1_18','CA20_18','CA21_18','CA22_18','CA23_18','CA24_18','CA2_18','CA3_18','CA4_18','CA5_18','CA6_18','CA7_18')
TA13 <- c('T10','T11','T13','T15','T16','T19','T2','T3','T5','T6','T7','T8')
TA18 <- c('TA11_18','TA12_18','TA13_18','TA14_18','TA15_18','TA16_18','TA17_18','TA18_18','TA19_18','TA1_18','TA20_18','TA21_18','TA22_18','TA23_18','TA24_18','TA2_18','TA3_18','TA4_18','TA5_18','TA6_18','TA7_18','TA8_18','TA9_18')
UL13 <- c('UL1','UL10','UL11','UL13','UL17','UL19','UL2','UL20','UL3','UL6','UL7','UL8','UL9')
UL18 <- c('UL10_18','UL11_18','UL12_18','UL13_18','UL14_18','UL15_18','UL16_18','UL17_18','UL18_18','UL19_18','UL1_18','UL20_18','UL21_18','UL22_18','UL23_18','UL24_18','UL2_18','UL3_18','UL4_18','UL5_18','UL6_18','UL7_18','UL8_18','UL9_18')
LL13 <- c('LL1','LL10','LL11','LL12','LL13','LL14','LL2','LL3','LL5','LL6','LL7','LL8','LL9')
LL18 <- c('LL12_18','LL13_18','LL14_18','LL15_18','LL16_18','LL17_18','LL18_18','LL19_18','LL1_18','LL20_18','LL21_18','LL22_18','LL23_18','LL24_18','LL2_18','LL3_18','LL4_18','LL5_18','LL6_18','LL7_18','LL8_18','LL9_18')
GH <- c('GH1','GH10','GH11','GH12','GH13','GH14','GH15','GH17','GH18','GH19','GH2','GH20','GH3','GH4','GH5','GH6','GH7','GH8','GH9')
vcf_file <- set.populations(vcf_file, list(GH,CA13,CA18,TA13,TA18,LL13,LL18,UL13,UL18), diploid=TRUE)
#vcf_file <- set.outgroup(vcf_file,new.outgroup=FALSE,diploid=TRUE)

### create windows ###
slide_vcf_file_75kb <- sliding.window.transform(vcf_file, 50000, 50000, type=2)

### do pop stats
slide_vcf_file_75kb <- F_ST.stats(slide_vcf_file_75kb, mode="nucleotide")
slide_vcf_file_75kb <- neutrality.stats(slide_vcf_file_75kb, FAST=TRUE)
slide_vcf_file_75kb <- diversity.stats.between(slide_vcf_file_75kb, nucleotide.mode=TRUE)
slide_vcf_file_75kb <- detail.stats(slide_vcf_file_75kb,site.spectrum=TRUE)


 calc_FST<-function(fst_list){
         tmp<-data.frame(GHP_CA13= fst_list[,1],
 		GH_CA18= fst_list[,2],
 		GH_TA13= fst_list[,3],
 		GH_TA18= fst_list[,4],
 		GH_LL13= fst_list[,5],
 		GH_LL18= fst_list[,6],
 		GH_UL13= fst_list[,7],
  		GH_UL18=fst_list[,8], 
  		CA13_CA18=fst_list[,9],
  		CA13_TA13=fst_list[,10],
 		CA13_TA18= fst_list[,11],
 		CA13_LL13= fst_list[,12],
 		CA13_LL18= fst_list[,13],
 		CA13_UL13= fst_list[,14],
 		CA13_UL18= fst_list[,15],
 		CA18_TA13= fst_list[,16],
  	CA18_TA18=fst_list[,17],
 		CA18_LL13=fst_list[,18],
 		CA18_LL18=fst_list[,19],
  	CA18_UL13= fst_list[,20],
 		CA18_UL18= fst_list[,21],
		TA13_TA18= fst_list[,22],
		TA13_LL13= fst_list[,23],
		TA13_LL18= fst_list[,24],
		TA13_UL13=fst_list[,25],
		TA13_UL18=fst_list[,26],
		TA18_LL13=fst_list[,27],
		TA18_LL18= fst_list[,28],
		TA18_UL13= fst_list[,29],
		TA18_UL18= fst_list[,30],
		LL13_LL18= fst_list[,31],
		LL13_UL13=fst_list[,32],
		LL13_UL18=fst_list[,33],
		LL18_UL13=fst_list[,34],
  	LL18_UL18=fst_list[,35],
  	UL13_UL18=fst_list[,36])

      if(length(tmp) == 0){
         return('NA')
         } else {
         tmp['chrom']<-rep(tid,nrow(tmp))
         tmp['window']<-1:nrow(tmp)
         tmp['window_start']<- (tmp$window-1)*50000
         tmp['window_end']<- tmp$window*50000
         tmp<-tmp[,c(37:40,1:36)]
         }
         }
FST_mat<-data.frame(calc_FST(t(slide_vcf_file_75kb@nuc.F_ST.pairwise)))
outfile_fst <-paste0("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/popgenome/",tid,".50kb.fst.popgenome.out")
 write.table(FST_mat, outfile_fst, sep = "\t", row.names=F, quote = F)

 calc_PI<-function(pi_list){
         tmp<-data.frame(GH = pi_list[,1]/50000,
                 CA13 = pi_list[,2]/50000,
                 CA18 = pi_list[,3]/50000,
                 TA13 = pi_list[,4]/50000,
		 TA18 = pi_list[,5]/50000,
		 LL13 = pi_list[,6]/50000,
		 LL18 = pi_list[,7]/50000,
		 UL13 = pi_list[,8]/50000,
		 UL18 = pi_list[,9]/50000)
         if(length(tmp) == 0){
         return('NA')
         } else {
         tmp['chrom']<-rep(tid,nrow(tmp))
         tmp['window']<-1:nrow(tmp)
         tmp['window_start']<- (tmp$window-1)*50000
         tmp['window_end']<- tmp$window*50000
         tmp<-tmp[,c(10:13,1:9)]
         }
 }
 
 PI_mat<-data.frame(calc_PI(slide_vcf_file_75kb@nuc.diversity.within))
 outfile_pi <-paste0("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/popgenome/",tid,".50kb.pi.popgenome.out")
 write.table(PI_mat, outfile_pi, sep = "\t", row.names=F, quote = F)

 calc_TD<-function(td_list){

         tmp<-data.frame(GH= td_list[,1],
                 CA13 = td_list[,2],
                 CA18 = td_list[,3],
  TA13 = td_list[,4],
		 TA18 = td_list[,5],
                 LL13 = td_list[,6],
                 LL18 = td_list[,7],
                 UL13 = td_list[,8],
                 UL18 = td_list[,9])
         if(length(tmp) == 0){
         return('NA')
         } else {
         tmp['chrom']<-rep(tid,nrow(tmp))
         tmp['window']<-1:nrow(tmp)
         tmp['window_start']<- (tmp$window-1)*50000
         tmp['window_end']<- tmp$window*50000
         tmp<-tmp[,c(10:13,1:9)]
         }
         }
 TD_mat<-data.frame(calc_TD(slide_vcf_file_75kb@Tajima.D))
 outfile_td <-paste0("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/popgenome/",tid,".50kb.td.popgenome.out")
 write.table(TD_mat, outfile_td, sep = "\t", row.names=F, quote = F)

}

x_vector<-seq(1,nrow(chr_length))
lapply(x_vector,chr_function)
