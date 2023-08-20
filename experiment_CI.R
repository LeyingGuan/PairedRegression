library(PREGS)
library(argparse)
source("simulations_helpers.R")
n=100
parser <- ArgumentParser(description='Simulation settings for PREGS')
parser$add_argument("--p", type = "integer", help = "num of features in Z")
parser$add_argument("--S", type = "numeric", help = "signal/power [0,.99]")
parser$add_argument("--D", type = "character", help = "design (AnovaBalance, Gaussian, Cauchy, CorGaussian,  AnovaUnbalance)")
parser$add_argument("--E", type = "character", help = "noise/error (gaussian, cauchy, exponential,  multinomial)")
args <- parser$parse_args()
args$D ="Gaussian"; args$E ="gaussian"
p = args$p
sig = args$S
cat( args$D, "\t", args$E, "\t", p, "\t", sig)

result_file_name = paste("CIevaluation",args$D, args$E, args$n, args$p, "signal", args$S*100,sep="_")
result_file_name = paste0(result_file_name, ".Rdata")
path = "/home/lg689/project/project/PairedRegression/Results/"


iter = 10; B = 2000; M = 2000
set.seed(2024)
seeds = sample(1:100000, iter, replace = F)
results = list()
results[["coverage"]] = data.frame(matrix(NA, ncol = 3), nrow = iter)
results[["length"]] =data.frame(matrix(NA, ncol = 3), nrow = iter)
colnames(results[["coverage"]]) <- colnames(results[["length"]]) = c("Normal", "Bootstrap", "Inversion")

alpha = 0.05

for(it in 1:iter){
  print(paste0("#########", it, "############"))
  dat = data_gen(n=n, p = p, M = M, design =args$D, noise = args$E,seed = seeds[it])
  if(as.numeric(sig) <=0.0){
    dat$beta = 0.0
  }else{
    dat$beta  = beta_gen(dat, power = sig)
  }
  dat$y = y_gen(dat$x, dat$beta, dat$epsMat)
  out_PREG <- PREGtest(x= dat$x, y = dat$y, z=dat$Z, B = B)
  CItable_list = PREGS_CI_table(out_PREG$res$Cmat)
  CI_PREG_invert = PREG_CI_construct(CItable_list, alpha = alpha)
  results[["coverage"]][it,3] = mean((CI_PREG_invert$beta_min <= dat$beta) & (CI_PREG_invert$beta_max >= dat$beta))
  results[["length"]][it,3] = median(CI_PREG_invert[,2]-CI_PREG_invert[,1])
  
  ###Normal approximation
  lm_fit = lm(dat$y~dat$x+dat$Z)
  lm_fit_summary = summary(lm_fit)
  qt_val <- qt(1 - alpha/2, df.residual(lm_fit))
  coef_se <- sapply(lm_fit_summary, function(z) z$coefficients[2,2])
  coef_est <- sapply(lm_fit_summary, function(z) z$coefficients[2,1])
  CI_lm = cbind(coef_est-coef_se*qt_val, coef_est+coef_se*qt_val)
  
  results[["coverage"]][it,1] = mean((CI_lm[,1] <= dat$beta) & (CI_lm[,2] >= dat$beta))
  results[["length"]][it,1] = median(CI_lm[,2]-CI_lm[,1])
  
  
  ###Bootstrap
  # CI_PREG_boot0 <- boot(dat, regress_func, R = 1000)
  # CI_PREG_boot = matrix(NA, nrow = ncol(dat$y), ncol = 2)
  # for(i in 1:ncol(dat$y)){
  #   regress_func <- function(dat, indices) {
  #     model <- lm(dat$y[indices,,drop = F] ~ dat$x[indices]+dat$Z[indices,,drop = F])
  #     return(coef(model)[2,])
  #   }
  #
  #   tmp = boot.ci(CI_PREG_boot0, index = i, type = "bca", conf =  1.0 -alpha)
  #   CI_PREG_boot[i,1]=tmp$bca[4]
  #   CI_PREG_boot[i,2]=tmp$bca[5]
  # }
  # mean((CI_PREG_boot[,1] <= dat$beta) & (CI_PREG_boot[,2] >= dat$beta))
  # median(CI_PREG_boot[,2]-CI_PREG_boot[,1])
  
  CI_PREG_boot = matrix(NA, nrow = ncol(dat$y), ncol = 2)
  ##get sufficient statistics
  
  for(i in 1:ncol(dat$y)){
    Xy <- cbind(cbind(dat$y[,i], dat$x),dat$Z)
    rfun <- function(Xy) {
      y <- Xy[, 1]
      model = lm(y~Xy[,-1])
      coef(model)[2]
    }
    tmp = bcajack(x = Xy, B = 500, func = rfun, m = 20, verbose = FALSE)
    CI_PREG_boot[i,1]=tmp$lims[1,1]
    CI_PREG_boot[i,2]=tmp$lims[9,1]
  }
  
  results[["coverage"]][it,2] =mean((CI_PREG_boot[,1] <= dat$beta) & (CI_PREG_boot[,2] >= dat$beta))
  results[["length"]][it,2] = median(CI_PREG_boot[,2]-CI_PREG_boot[,1])
  
  print(results[["coverage"]][it,])
  
  saveRDS(results, file = paste0(path,result_file_name))
}
