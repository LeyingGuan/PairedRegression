#!/bin/bash

#SBATCH --job-name=example_job
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL

module load miniconda/23.3.1
conda init bash
conda activate /home/lg689/project/project/conda_env/PREGSenv
Rscript experiments_sim.R --n 100 --p 2 --D AnovaBalance --E gaussian --S 0.0

