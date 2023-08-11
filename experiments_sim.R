#design:"AnovaBalance", "Gaussian", "Cauchy", "CorGaussian",  "AnovaUnbalance"
#noise: "gaussian", "cauchy", "exponential", "multinomial"
#beta = 0; beta = 80% power for Gaussian
#alpha = 0.05, 0.01, 0.001
# n = 100; ncol(Z)=2, 5, 10 (intercept included)
#Each run: B = 2000, M = 2000, repeat 100 times
# Methods for comparison: PRegs, FL, Permutation, t-test, cyclic
library(PREGS)
library(argparse)
source("simulations_helpers.R")


messages=c("sample size", "num of features in Z", "design", "noise/error", "signal/power",
           "run CPT")
parser <- ArgumentParser(description='Simulation settings for PREGS')
parser$add_argument("--n", type = "integer", help = "sample size")
parser$add_argument("--p", type = "integer", help = "num of features in Z")
parser$add_argument("--D", type = "character", help = "design (AnovaBalance, Gaussian, Cauchy, CorGaussian,  AnovaUnbalance)")
parser$add_argument("--E", type = "character", help = "noise/error (gaussian, cauchy, exponential,  multinomial)")
parser$add_argument("--S", type = "numeric", help = "signal/power [0,.99]")
parser$add_argument("--CPT", type = "logical", help = "If run CPT")
args <- parser$parse_args()
#args$D = "Cauchy"; args$E = "multinomial"; args$S=0; args$p=20;args$n=100

cat(messages[1], args$n, "\n")
cat(messages[2], args$p, "\n")
cat(messages[3], args$D, "\n")
cat(messages[4], args$E, "\n")
cat(messages[5], args$S, "\n")
cat(messages[6], args$CPT, "\n")


iter = 20; B = 2000; M = 2000
set.seed(2024)
seeds = sample(1:100000, iter, replace = F)
results = list()
results[["pvalues"]] = list()
results[["setting"]] =args
if(args$S<=0.0){
  result_file_name = paste(args$D, args$E, args$n, args$p, "null", sep="_")
}else{
  result_file_name = paste(args$D, args$E, args$n, args$p, "alternative",sep="_")
}
if(args$CPT){
  result_file_name=paste0(result_file_name, "_", "wiCPT",".Rdata")
}else{
  result_file_name=paste0(result_file_name, "_", "woCPT",".Rdata")
}
path = "/home/lg689/project/project/PairedRegression/Results/"

for(it in 1:iter){
  print(paste0("#########", it, "############"))
  dat = data_gen(n=as.integer(args$n), p = as.integer(args$p), M = M, design =args$D, noise = args$E,seed = seeds[it])
  if(as.numeric(args$S) <=0.0){
    dat$beta = 0.0
  }else{
    dat$beta  = beta_gen(dat$x, dat$Z, dat$epsMat, power = as.numeric(args$S))
  }
  dat$y = y_gen(dat$x, dat$beta, dat$epsMat)
  results[["pvalues"]][[it]] = sim_comparisons_singleSetting(dat, B = B, run_CPT=run_CPT)

  print(apply(results[["pvalues"]][[it]]<=0.05,2,mean))
  #print(apply(results[["pvalues"]][[it]]<=0.05,2,mean))
  #print(apply(results[["pvalues"]][[it]]<=0.05,2,mean))
  saveRDS(results, file = paste0(path,result_file_name))
}

