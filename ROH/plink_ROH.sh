#!/bin/bash
#SBATCH -D .
#SBATCH -p pq
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=test
#SBATCH --error=test.err.txt
#SBATCH --output=test.out.txt
#SBATCH --export=All


module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0


maf_vcf=/gpfs/ts0/projects/Research_Project-T109423/ali/FIBR_transfer/ftp1.sequencing.exeter.ac.uk/V0085/11_trimmed/PopGenome/FIBR_2018_liftover_phased.vcf.gz
pop_map_GH=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/GH_popmap.txt
pop_map_CA13=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/CA_2013_popmap.txt
pop_map_CA18=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/CA_2018_popmap.txt
pop_map_TA13=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/TA_2013_popmap.txt
pop_map_TA18=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/TA_2018_popmap.txt
pop_map_UL13=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/UL_2013_popmap.txt
pop_map_UL18=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/UL_2018_popmap.txt
pop_map_LL13=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/LL_2013_popmap.txt
pop_map_LL18=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/vcftools/popmaps/LL_2018_popmap.txt

#subset for each POP
#then run plink on each
vcftools --gzvcf $maf_vcf --keep $pop_map_GH  --recode --remove-filtered-all --out GH_only

/lustre/bin/plink --vcf $DIR/GH_phased_autosomes.recode.vcf.gz --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out $OUT/GHP_500K_1het --allow-extra-chr

~/plink --vcf GH_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out GHP_500K_1het --allow-extra-chr


vcftools --gzvcf $maf_vcf --keep $pop_map_CA13  --recode --remove-filtered-all --out CA13_only

~/plink --vcf CA13_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out CA13_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_CA18  --recode --remove-filtered-all --out CA18_only

~/plink --vcf CA18_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out CA18_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_TA13  --recode --remove-filtered-all --out TA13_only

~/plink --vcf TA13_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out TA13_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_TA18  --recode --remove-filtered-all --out TA18_only

~/plink --vcf TA18_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out TA18_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_LL13  --recode --remove-filtered-all --out LL13_only

~/plink --vcf LL13_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out LL13_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_LL18  --recode --remove-filtered-all --out LL18_only

~/plink --vcf LL18_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out LL18_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_UL13  --recode --remove-filtered-all --out UL13_only

~/plink --vcf UL13_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out UL13_500K_1het --allow-extra-chr

vcftools --gzvcf $maf_vcf --keep $pop_map_UL18  --recode --remove-filtered-all --out UL18_only

~/plink --vcf UL18_only.recode.vcf --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out UL18_500K_1het --allow-extra-chr

