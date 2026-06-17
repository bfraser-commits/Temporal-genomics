#!/bin/bash
#SBATCH -p pq
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -A Research_Project-T110748
#SBATCH --job-name=fibr_slim
#SBATCH --error=fibr_slim.err.txt
#SBATCH --output=fibr_slim.out.txt
#SBATCH --export=All
#SBATCH --array=1-100%25

module load GCCcore/13.2.0

/lustre/home/bf299/build/slim  -d OUTPUT_PREFIX="'vcf_outs_mult_pops/run_${SLURM_ARRAY_TASK_ID}'" input_vcf_output_vcf_mult_outs_all_pops_FINAL.slim
