#!/bin/bash

# Define the set of character values
declare -a design_values=("AnovaBalance" "Gaussian" "Cauchy" "Paired")
declare -a noise_values=("gaussian" "cauchy")
declare -a power_values=(0.0 0.3 0.5 0.7 0.9)
declare -a n_values=(100)
declare -a p_values=(2 6 16)
FILE=$1
for design in "${design_values[@]}"
do
  for noise in "${noise_values[@]}"
  do
    for power in "${power_values[@]}"
    do
      for p in "${p_values[@]}"
      do
        echo "module load miniconda/23.3.1;conda activate /home/lg689/project/project/conda_env/PREGSenv;Rscript experiment_CI.R --p $p --D $design --E $noise --S $power" >> $FILE
      done
    done
  done
done
