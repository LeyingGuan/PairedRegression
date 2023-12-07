library(PREGS)
library(MASS)
source("dataGen.R")
source("CPT_RPT.R")

#source("/Users/lg689/Downloads/CPT-master/code/expr_functions.R")
sim_comparisons_singleSetting = function(dat, run_CPT=TRUE, B = 2000){
  
  pvals_collection = data.frame(matrix(NA, ncol = 5, nrow = ncol(dat$y)))
  colnames(pvals_collection) = c("Ftest", "PERMtest", 
                                 "FLtest", "PAMLRT","CPT")
  ###t test
  fitted_ttest = summary(lm(dat$y~dat$x+dat$Z))
  pvals_collection[,1] = as.vector(sapply(fitted_ttest, function(z) z$coefficients[2,4]))
  
  ##simple permutations
  pvals_collection[,2] = PERMtest(x= dat$x, y = dat$y, z=dat$Z, B = B)$res$unsigned
  ###FL-test (F)
  pvals_collection[,3] =  FLtest(x= dat$x, y = dat$y, z=dat$Z, B = B)$res$unsigned
  
  ## PRegs
  pvals_collection[,4] = PREGtest(x= dat$x, y = dat$y, z=dat$Z, B = B)$res$unsigned


  ##cyclic permutations
  ### genetic alg ordering
  if( run_CPT){
    rounds <- 100 + 900; m = 19;M = NULL
    X0 = cbind(dat$x, dat$Z[,-ncol(dat$Z)])
    X = cbind(dat$x, dat$Z)
    res1 <- find_eta_GA(X0, m, testinds = 1, popSize = 10, rounds = rounds,M = M)
    ordering = res1$ordering
    for(l in 1:ncol(dat$y)){
      out0 = CPT(dat$y[,l], X, ordering = ordering,testinds = 1,alpha = 0.05,returnCI = FALSE)
      pvals_collection[l,5]=  out0$pval[1]
    }
  }
  # apply(pvals_collection<=0.05,2,mean)
  # apply(pvals_collection<=0.01,2,mean)
  # apply(pvals_collection<=0.001,2,mean)

  return(pvals_collection)
  
}
