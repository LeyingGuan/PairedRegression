library(PREGS)
library(MASS)
library(plyr)


generatePermutations <- function(n, K){
  ind <- head(sample(n), n - n %% (K+1))
  mx <- matrix(ind, nrow=K+1)
  perm_ind <- matrix(0, K, n)
  for (r in 1:K){
    perm_ind[r, ] <- plyr::mapvalues(1:n, mx, mx[c((r+1):(K+1), 1:r), ])
  }
  return(perm_ind)
}

# PRT (code modified from Wang et al 2023)
RPT <- function(X, Y, Z, K=99){
  n <- dim(X)[1]; p <- dim(X)[2]
  if(!is.null(dim(Y))){
    stat1 <- stat2 <- matrix(0,ncol=K, nrow = ncol(Y))
      rep(0, K)
  }else{
    stat1 <- stat2 <- rep(0, K)
  }

  perm_ind <- generatePermutations(n, K)

  for (r in 1:K){
    idx <- perm_ind[r, ]  # P is a perm matrix with P[i, idx[i]] = 1
    VtildeVtilde <- tryCatch({tmp <- cbind(X, X[idx,]); diag(n) - tmp %*% solve(t(tmp) %*% tmp, t(tmp))}, error=function(e){print("intertible error for tmp"); "error"})
    if(identical(VtildeVtilde, "error")){
      #there may be overlap between spaces of X and X[idx, ] or X may not be full rank
      tmp <- svd(cbind(X, X[idx,]), nu=2*p)$u
      VtildeVtilde <- diag(n) - tmp %*% t(tmp)
    }
    #vectors <- svd(tmp, nu=n)$u[, -(1:(2 * p))]
    Vz <- as.vector(VtildeVtilde %*% Z)
    if(is.null(dim(Y))){
      stat1[r] <- as.numeric(sum(Vz * Y))
      stat2[r] <- as.numeric(sum(Vz * Y[idx]))
    }else{
      stat1[,r] <- as.numeric(apply(Y, 2, function(z){sum(Vz*z)}))
      stat2[,r] <- as.numeric(apply(Y, 2, function(z){sum(Vz*z[idx])}))
    }
  }
  if(is.null(dim(stat1))){
    stat1 = matrix(stat1,nrow = 1)
    stat2 = matrix(stat2,nrow = 1)
  }
  pval <-sapply(1:nrow(stat1), function(i){ (sum(abs(stat2[i,]) >= abs(stat1[i,])) + 1) / (K+1)})
  pval_theory <-sapply(1:nrow(stat1), function(i){ (sum(abs(stat2[i,]) >= min(abs(stat1[i,]))) + 1) / (K+1)})
  return(list(pval=pval, pval_theory=pval_theory))
}

data_gen = function(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1, m =1){
  xZ = matrix(NA, ncol = p+1, nrow = n)
  if(design == "AnovaBalance"){
    #xZ = t(rmultinom(n = n, size = 1, prob = rep(1.0/(p+1), p+1)))
    m = floor(n/(p+1))
    xZ = diag(1, p+1)
    #previously rbind(xZ, m)
    for(l in 1:(m-1)){
      xZ = rbind(xZ, diag(1, p+1))
    }
    if(nrow(xZ)<n){
      xZ = rbind(xZ, matrix(0, ncol = p+1, nrow = n-nrow(xZ)))
    }
    xZ[,p+1]=1
  }else if(design == "Gaussian"){
    xZ = matrix(rnorm(n*(p+1)), ncol = p+1)
    
  }else if(design == "Cauchy"){
    xZ = matrix(rcauchy(n*(p+1)), ncol = p+1)
    
  }else if(design == "CorGaussian"){
    xZ[,-1] = matrix(rnorm(n*p), ncol = p)
    xZ[,p+1] = 0
    ww = rnorm(p-1)
    ww = ww/sqrt(sum(ww**2))
    xZ[,1] = rnorm(n)
    for(l in 1:length(ww)){
      xZ[,1] = xZ[,1]+ww[l]*xZ[,1+l]
    }
    xZ[,1+p]=1.0
  }else if(design == "AnovaUnbalance"){
    xZ[,]=0
    for(j in 1:p){
      xZ[j,j]=1
    }
    xZ[,p+1]=1
  }else if(design == "Paired"){
    xZ = diag(1, p+1)
    #previously rbind(xZ, m)
    for(l in 1:m){
      xZ = rbind(xZ, rep(1, p+1))
    }
    if(nrow(xZ)<n){
      xZ = rbind(xZ, matrix(0, ncol = p+1, nrow = n-nrow(xZ)))
    }
    xZ[,p+1]=1
  }
  else{
    stop("unknown design")
  }
  xZ = apply(xZ, 2, function(z) z/sqrt(sum(z*z)))
  x = xZ[,1]
  Z = xZ[,-1]
  Z[,p] = 1
  epsMat = matrix(NA, ncol = M, nrow = n)
  if(noise == "gaussian"){
    epsMat = matrix(rnorm(n*M), ncol = M)
  }else if(noise == "cauchy"){
    epsMat = matrix(rcauchy(n*M), ncol = M)
  }else if(noise == "exponential"){
    epsMat = matrix(rexp(n*M)-1, ncol = M)
  }else if(noise == "multinomial"){
    epsMat = matrix(rnorm(n*M), ncol = M)
    for(l in 1:ncol(epsMat)){
      idx = sample(1:n,1, replace = F)
      epsMat[idx,l]= epsMat[idx,l]+2*(rbinom(length(idx), size = 1, prob = .5)-.5)*max(n*n, 1E4)
    }
  }else{
    stop("unknown noise")
  }
  x = matrix(x, ncol = 1)
  return(list(x = x, Z = Z, epsMat = epsMat))
}

beta_gen = function(dat, beta_grids =2^(c(1:400)/30),power = .9){
  x = dat$x
  Z = dat$Z
  epsMat = dat$epsMat
  n = nrow(Z)
  fitted =lm(epsMat ~ x+Z) 
  rx = lm(x~Z)$residuals
  reps_ = fitted$residuals
  reps = lm(epsMat~Z)$residuals
  sigma_mat_resid = rx%o%beta_grids
  residuals_denominator =apply(reps_^2, 2, sum)/(n-1-ncol(Z))
  #resid1 = reps[,1] +sigma_mat_resid[,which(beta_grids==beta)]
  ##y2=epsMat[,1]+x*beta
  #resid2 =  lm(y2~Z)$residuals
  #all.equal(resid1, resid2)
  residuals_numerator = sapply(1:length(beta_grids), function(i){
    apply((reps +sigma_mat_resid[,i])^2-reps_^2,2,sum)
  })
  Fstats = residuals_numerator/residuals_denominator
  pvals = pf(Fstats, df1 = 1, df2 = (n-1-ncol(Z)),lower.tail = F)
  powers = apply(pvals<=0.05,2,mean)
  id = which.min(abs(powers-power))
  # names(powers) = as.character(beta_grids)
  # print(powers)
  beta = beta_grids[id]
  return(beta)
}

y_gen = function(x, beta, epsMat){
  y = apply(epsMat, 2, function(z) z+x*beta)
  return(y)
}


# dat = data_gen(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1)
# dat$beta  = beta_gen(dat$x, dat$Z, dat$epsMat, power = .9)
# dat$y = y_gen(dat$x, dat$beta, dat$epsMat)

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
