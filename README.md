# Permutation-Augmented Linear Model Regression Test (PAMLRT)
A conformal test of non-zero coefficient in linear models via permutation-augmentation.

[https://arxiv.org/pdf/2309.05482.pdf]



## Package installation

[Note that the package name is PREGS instead of PAMLRT for now]

library(devtools)

install_github("LeyingGuan/PairedRegression/PREGS")

## Example



## Manuscript reproducibility

## Comparisons of coverages and powers using different methods

experiments_sim.R: The R script experiments_sim.R runs PALMRT (PREGS), CPT, RPT, Ftest, FLtest, PERMtest for a given design, noise distribution, feature dimension, signal strength.

simulaiton_helpers.R: Helper functions for data generation ( design, noise distribution, feature dimension, signal strength), runing experiment in experiments_sim.R.

CPT_RPT.R: codes for running CPT and RPT (modified code from the original CPT and RPT articles).

## Confidence interval comparisons.

