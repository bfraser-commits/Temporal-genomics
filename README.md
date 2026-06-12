# Temporal genomics code
#code to accompany FIBR temporal paper

#Read mapping, SNP calling, and filtering were done as in Whiting et al. 2020

https://github.com/josieparis/gatk-snp-calling

#Demographic changes across time figure

census_plot_new.R

#PCA analysis

pca_FIBR_temp_plink.sh

#figure was plotted using R

pca_plot.R

#Whole genome population genomics

#popgenome script:

popgenome_FIBR_temp_50kb.R

#R plotting and analysis scripts:

FST_year_figure.R

pi_figures.R

TajD_figures.R

#ROH analysis

#plink script

plink_ROH.sh

#R plotting and analysis script

ROH_analysis.R

#SLIMs for neutral analysis

#SLIM analysis

input_vcf_output_vcf_mult_outs_all_pops_FINAL.slim

#run multiple times

run_vcf_all_pops.sh 

#CVTK analysis

#first make a vcf that has GH replicated for the design in cvtk to work:

#cvtk modified python script:

FIBR_real_data.py

#to run on each chromosome

run_cvtk_chrs.sh

FIBR_chr_test.py

chr_plots_violin.R

#to run cvtk on slims

run_cvtk_slims.sh 

slims_out_all.sh 

#to plot results

plot_WG_percentile_calc.R

diagnostic_plots.R

#Analyse and find cut-offs for windows

#will take input and windowise it for comparisons

windowise_slims.R

#will take slim output to get AF cut-offs

explore_data_nopol_13_18_$POP.R

#will compare observed windows to simulate ones, get cutoffs and do overlap analysis and sup figures

window_analysis_slim_vcf_cutoff.R

#to plot and analyse AF direction of change 

#script to pull out annotations and run through GO analysis

candidate_region_analysis.sh

Biomart.R 




