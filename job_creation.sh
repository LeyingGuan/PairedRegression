#!/bin/bash

# Define the set of character values
declare -a design_values=("AnovaBalance" "Gaussian" "Cauchy" "AnovaUnbalance")
declare -a noise_values=("gaussian" "cauchy" "exponential" "multinomial")
declare -a power_values=(0.0 0.9)
declare -a n_values=(100 500)
declare -a p_values=(2 5 20)
FILE=$1
for design in "${design_values[@]}"
do
  for noise in "${noise_values[@]}"
  do
    for power in "${power_values[@]}"
    do
      for p in "${p_values[@]}"
      do
        for n in "${n_values[@]}"
	do
          echo "module load miniconda/23.3.1;conda activate /home/lg689/project/project/conda_env/PREGSenv;Rscript experiments_sim.R --n 100 --p $p --D $design --E $noise --S $power" >> $FILE
        done
      done
    done
  done
done
