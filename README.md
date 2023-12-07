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

CPT_RPT.R: codes for running CPT and RPT (modified code from the original CPT and RPT articles).

simulaiton_helpers.R: Helper functions for data generation ( design, noise distribution, feature dimension, signal strength), runing experiment in experiments_sim.R.

experiments_sim.R: The R script experiments_sim.R runs PALMRT (PREGS), CPT, RPT, Ftest, FLtest, PERMtest for a given design, noise distribution, feature dimension, signal strength. 

Below are the command lines used to generate all results related to the covarage and power comparison experiments.

```ruby
' Data generation

' Bulk job array creation and dSQ job running

' Example of running experiments_sim.R for one setting
```




## Confidence interval comparisons.

