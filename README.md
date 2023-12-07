# Permutation-Augmented Linear Model Regression Test (PAMLRT)
A conformal test of non-zero coefficient in linear models via permutation-augmentation.

[https://arxiv.org/pdf/2309.05482.pdf]



## Package installation

[Note that the package name is PREGS instead of PAMLRT for now]

library(devtools)

install_github("LeyingGuan/PairedRegression/PREGS")

## Example use PALMRT (PREGS).
dataGen.R: generate synthetic data.

```ruby
## M: dimension of y ( M independent realization of n noise)
dat = data_gen(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1)

## add signal to y
dat$beta  = beta_gen(dat$x, dat$Z, dat$epsMat, power = .9)
dat$y = y_gen(dat$x, dat$beta, dat$epsMat)

B = 2000
## x: n by 1; y: n by M; Z: n by p (note that one of the p column is the vector of 1's, so there are actually (p -1) covariates to be adjusted for).
PREGout = PREGtest(x= dat$x, y = dat$y, z=dat$Z, B = B)

## pvalues for unsigned test (F-test), signed version based on coefficients is also available
pvals = PREGout$res$unsigned

## Construct CI based on inversion

CItable_list = PREGS_CI_table(PREGout$res$Cmat)

CI_PREG_invert = PREG_CI_construct(CItable_list, alpha = alpha) 


```

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

