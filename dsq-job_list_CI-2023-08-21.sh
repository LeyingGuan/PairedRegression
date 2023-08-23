#!/bin/bash
#SBATCH --output dsq-job_list_CI-%A_%3a-%N.out
#SBATCH --array 0-119
#SBATCH --job-name dsq-job_list_CI
#SBATCH --mem-per-cpu 4g -t 8:00:00 --mail-type ALL

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/project/guan_leying/lg689/project/PairedRegression/PairedRegression/job_list_CI.txt --status-dir /gpfs/gibbs/project/guan_leying/lg689/project/PairedRegression/PairedRegression

