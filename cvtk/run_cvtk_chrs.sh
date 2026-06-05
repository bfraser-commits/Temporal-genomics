#!/bin/bash

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
FAI=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/STAR.holi11.analysis.genome.fasta.fai

#for CHR in $(cat $FAI | awk '$2 > 1000000' | cut -f1)
for CHR in chr21 chr22 chr23

do

  # Construct the paths
  vcf="/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/cvtk/merged_output_unphase.vcf.gz"
  popmap="/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/cvtk/real_data/merged_popmap_int.txt"	
  output="${CHR}_out"

  vcftools --gzvcf $vcf --chr "${CHR}" --recode --remove-filtered-all --out ${CHR}_
  # Run the Python script with the constructed paths
  python FIBR_chr_test.py --vcf ${CHR}_.recode.vcf --popmap "$popmap" --output "$output"
done
