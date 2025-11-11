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
module load PLINK/2.00-alpha1-x86_64

maf_vcf=/gpfs/ts0/projects/Research_Project-T109423/ali/FIBR_transfer/ftp1.sequencing.exeter.ac.uk/V0085/11_trimmed/PopGenome/FIBR_2018_liftover_phased.vcf.gz


#run plink
plink2 --vcf $maf_vcf --pca --indep-pairwise 50 5 0.2 --allow-extra-chr --out FIBR_thin_plink
