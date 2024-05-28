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
require(PREGS)
source("dataGen.R")
dat = data_gen(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1)

## add signal to y
dat$beta  = beta_gen(dat, power = .9)
dat$y = y_gen(dat$x, dat$beta, dat$epsMat)

B = 2000
## x: n by 1; y: n by M; Z: n by p (note that one of the p column is the vector of 1's, so there are actually (p -1) covariates to be adjusted for).
PREGout = PREGtest(x= dat$x, y = dat$y, z=dat$Z, B = B)

## pvalues for unsigned test (F-test), signed version based on coefficients is also available, but the CI construction is not yet implemented
pvals = PREGout$res$unsigned

## Construct CI based on inversion
### Cmat is the B by 4 by M Cmatrix in the CI section [PALMRT paper], used for constructing CI for each reliazation of y[,l], l=1,...,M
### CItable_list is a length M list, with element in each list correponds to  CI-level lookup table for one (x, Z, y[,l]). The table contains many intermediate quantities as described in the paper.

CItable_list = PREGS_CI_table(PREGout$res$Cmat)

###construct CI for a given level; CI_PREG_invert is M by 2 (lower boundary, higher boundary)
alpha = 0.05
CI_PREG_invert = PREG_CI_construct(CItable_list, alpha = alpha) 

```

## Manuscript reproducibility

## Comparisons of coverages and powers using different methods

CPT_RPT.R: codes for running CPT and RPT (modified code from the original CPT and RPT articles).

simulaiton_helpers.R: Helper functions for data generation ( design, noise distribution, feature dimension, signal strength), runing experiment in experiments_sim.R.

experiments_sim.R: The R script experiments_sim.R runs PALMRT (PREGS), CPT, RPT, Ftest, FLtest, PERMtest for a given design, noise distribution, feature dimension, signal strength. 

Note that CPT and PRT algorithms require some other additional packages (for CPT, R package gaoptim for genetic algorithm is no longer available on R CRAN for R >=4.2.3, and manual installation is needed.) 

### Run experiments_sim.R for one setting
```ruby
## --n=sample size; --p=ncol(Z)+1 (suppose Z does not include intercept); --D=design; --E=noise distribution; --S=signal strength (measured by power), 0=null, 1=100% power; --CPT=if run CPT
Rscript experiments_sim.R --n 100 --p 2 --D AnovaBalance --E gaussian --S 0.0 --CPT TRUE
```
### Run all analysis
Below are the command lines used to generate all results related to the covarage and power comparison experiments.
```ruby
## create job list
bash job_creation_wiCPT.sh experiments_main_job.txt

## large-scale experiments with dSQ
dsq --job-file experiments_main_job.txt [....other arguments]
dsq-experiments_main_job-2023-12-06.sh
```
### Figure Creation

## Experiments for CI quality evaluation




## Confidence interval comparisons.

