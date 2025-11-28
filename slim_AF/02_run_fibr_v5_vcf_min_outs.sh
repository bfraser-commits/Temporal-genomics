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

# Set environment
#module load GCC/8.1.0-2.30
module load GCCcore/13.2.0
#module load Perl/5.26.1-foss-2018a

#MAIN=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/slims_2018
#cd $MAIN

# Set up metadata
fibr_pops=(CA TY LL UL)


# Make somewhere to store results
NO_SELECTION=outputs_vcf_min_outs

# Run a batch array for all slim sims


/lustre/home/bf299/build/slim  -d 'demo_data="data/CA_2018_simulation_demography_min.txt"' input_vcf_output_pos_slim.slim > $NO_SELECTION/CA_test_run_iter${SLURM_ARRAY_TASK_ID}.res
/lustre/home/bf299/build/slim -d 'demo_data="data/TY_2018_simulation_demography_min.txt"' input_vcf_output_pos_slim.slim > $NO_SELECTION/TY_test_run_iter${SLURM_ARRAY_TASK_ID}.res
/lustre/home/bf299/build/slim -d 'demo_data="data/UL_2018_simulation_demography_min.txt"' input_vcf_output_pos_slim.slim > $NO_SELECTION/UL_test_run_iter${SLURM_ARRAY_TASK_ID}.res
/lustre/home/bf299/build/slim -d 'demo_data="data/LL_2018_simulation_demography_min.txt"' input_vcf_output_pos_slim.slim  > $NO_SELECTION/LL_test_run_iter${SLURM_ARRAY_TASK_ID}.res
