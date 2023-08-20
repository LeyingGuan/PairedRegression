#!/bin/bash
#SBATCH --output dsq-experiments_main_job-%A_%3a-%N.out
#SBATCH --array 0-599
#SBATCH --job-name dsq-experiments_main_job
#SBATCH --mem-per-cpu 1g -t 2:00:00 --mail-type ALL

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/project/guan_leying/lg689/project/PairedRegression/PairedRegression/experiments_main_job.txt --status-dir /gpfs/gibbs/project/guan_leying/lg689/project/PairedRegression/PairedRegression

