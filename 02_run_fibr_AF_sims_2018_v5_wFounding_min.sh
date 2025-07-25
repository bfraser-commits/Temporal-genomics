#!/bin/bash
#SBATCH -p pq
#SBATCH --time=02:00:00
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

MAIN=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/FIBR_temporal/slims_2018
#cd $MAIN

# Set up metadata
fibr_pops=(CA TY LL UL)
#fibr_start_sizes=(104 128 76 76)


# # Make somewhere to store results
# WITH_SELECTION=outputs/FIBR_runs_selection
# mkdir $WITH_SELECTION
# ALPHA=0.05
#
# # Run a batch array for all slim sims
# slim -d founding_size="${fibr_start_sizes[0]}" -d alpha=$ALPHA -d 'demo_data="data/CA_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction_with_selection.slim > $WITH_SELECTION/CA_test_run_iter${SLURM_TASK_ARRAY_ID}.res
# slim -d founding_size="${fibr_start_sizes[1]}" -d alpha=$ALPHA -d 'demo_data="data/TY_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction_with_selection.slim > $WITH_SELECTION/TY_test_run_iter${SLURM_TASK_ARRAY_ID}.res
# slim -d founding_size="${fibr_start_sizes[2]}" -d alpha=$ALPHA -d 'demo_data="data/LL_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction_with_selection.slim > $WITH_SELECTION/UL_test_run_iter${SLURM_TASK_ARRAY_ID}.res
# slim -d founding_size="${fibr_start_sizes[3]}" -d alpha=$ALPHA -d 'demo_data="data/UL_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction_with_selection.slim > $WITH_SELECTION/LL_test_run_iter${SLURM_TASK_ARRAY_ID}.res

# Make somewhere to store results
NO_SELECTION=outputs_test
#mkdir $NO_SELECTION

# Run a batch array for all slim sims
#NE is roughly 20,000 from Whiting et al

/lustre/home/bf299/build/slim -d founding_size="20000" -d 'demo_data="data/CA_2018_simulation_demography_min.txt"' -d 'burnin_path="GH_treeseq_out.trees"' slim/simulate_fibr_introduction_with_sampling_err_v5_wFounding.slim  > $NO_SELECTION/CA_test_run_iter${SLURM_ARRAY_TASK_ID}.res
/lustre/home/bf299/build/slim -d founding_size="20000" -d 'demo_data="data/TY_2018_simulation_demography_min.txt"' -d 'burnin_path="GH_treeseq_out.trees"' slim/simulate_fibr_introduction_with_sampling_err_v5_wFounding.slim  > $NO_SELECTION/TY_test_run_iter${SLURM_ARRAY_TASK_ID}.res
/lustre/home/bf299/build/slim -d founding_size="20000" -d 'demo_data="data/LL_2018_simulation_demography_min.txt"' -d 'burnin_path="GH_treeseq_out.trees"' slim/simulate_fibr_introduction_with_sampling_err_v5_wFounding.slim  > $NO_SELECTION/LL_test_run_iter${SLURM_ARRAY_TASK_ID}.res
/lustre/home/bf299/build/slim -d founding_size="20000" -d 'demo_data="data/UL_2018_simulation_demography_min.txt"' -d 'burnin_path="GH_treeseq_out.trees"' slim/simulate_fibr_introduction_with_sampling_err_v5_wFounding.slim  > $NO_SELECTION/UL_test_run_iter${SLURM_ARRAY_TASK_ID}.res
