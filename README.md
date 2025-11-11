# FIBR_temporal
#code to accompany FIBR temporal paper
#Read mapping, SNP calling, and filtering were done as in Whiting et al. 2020
#cite Paris github
#Demographic changes across time figure: R script

#PCA analysis
#plink scripts:
#figure was plotted using R

#Whole genome population genomics
#popgenome scripts
#R plotting script

#ROH analysis
#plink script
#R plotting script

#CVTK analysis
#cvtk modified python script
#R plotting script

#Simulation of Selective sweeps
#First set of slim analysis uses a forward WF model coalescent burn-in, a minimum population size per generation and does not filter for low frequencies in the output

02_run_fibr_AF_sims_2018_v5_wFounding_min.sh

#calls this script

slim/simulate_fibr_introduction_with_sampling_err_v5_wFounding.slim

#using a non-polarised non-MAF summary of the data

explore_data_nopol_13_18.R

#to make supplemental plots, cut-offs, and windowise

~/Documents/FIBR_temporal/slims_test/compare_slim_obs_straight_cutoff.R

#overlap of outlier windows

~/Documents/FIBR_temporal/compare_window_outliers.R

####Second set of slim analysis uses forward WF model with a vcf from GH of 1 chromosome, all other parameters are the same as above, this is likely what will be going in the paper in this directory: people/bonnie/FIBR_temporal/vcftools/slims/

02_run_fibr_v5_vcf_min_outs.sh

#calls:
input_vcf_output_pos_slim.slim

#will take input and windowise it for comparisons
windowise_slims.R

#will take slim output to get AF cut-offs
explore_data_nopol_13_18_$POP.R

#will compare observed windows to simulate ones, get cutoffs and do overlap analysis and sup figures
~/Documents/FIBR_temporal/window_analysis_slim_vcf_cutoff.R

#to plot and analyse AF direction of change 

#annotations of genome
#cite Paris github

#script to pull out annotations

#to do GO enrichment and gene annotation of candidates




