library(PREGS)
library(MASS)

data_gen = function(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1){
  xZ = matrix(NA, ncol = p+1, nrow = n)
  if(design == "AnovaBalance"){
    #xZ = t(rmultinom(n = n, size = 1, prob = rep(1.0/(p+1), p+1)))
    m = floor(n/(p+1))
    xZ = diag(1, p+1)
    for(l in 1:m){
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
  }else{
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

beta_gen = function(x, Z,epsMat, power = .9){
  n = nrow(Z)
  rx = lm(x~Z)$residuals
  reps_ = lm(epsMat ~ x+Z)$residuals
  reps = lm(epsMat~Z)$residuals
  s2x = sum(rx^2)
  sigma2= apply(reps_^2,2,sum)/(n-1-ncol(Z))
  #80% power for 0.05
  df = n-1-ncol(Z)
  q1= abs(qt(p = 0.025, df = df))
  q2A = qt(p=power,df = df,lower.tail = F)
  beta = mean(sigma2)*(q1-q2A)/sqrt(sum(rx^2))
  
  #q1= abs(qt(p = 0.001/2, df = df))
  #pnorm(q1 - beta/mean(sigma2), lower.tail = F)
  return(beta)
  
}

y_gen = function(x, beta, epsMat){
  y = apply(epsMat, 2, function(z) z+x[,1]*beta)
  return(y)
}


# dat = data_gen(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1)
# dat$beta  = beta_gen(dat$x, dat$Z, dat$epsMat, power = .9)
# dat$y = y_gen(dat$x, dat$beta, dat$epsMat)

#source("/Users/lg689/Downloads/CPT-master/code/expr_functions.R")
sim_comparisons_singleSetting = function(dat, B = 2000){
  ###t test
  fitted_ttest = summary(lm(dat$y~dat$x+dat$Z))
  pval_ttest = as.vector(sapply(fitted_ttest, function(z) z$coefficients[2,4]))
  
  ###FL-t test
  pval_FL_twosided = rep(NA, ncol(dat$y))
  pval_FL_onesided = rep(NA, ncol(dat$y))
  for(l in 1:ncol(dat$y)){
    y = dat$y[,l]
    x = dat$x
    z = dat$Z
    tmp1partial <- FL(y=y, x=x, z=z, add_intercept = T, B = B, statType = "partial")$res
    pval_FL_twosided[l] = tmp1partial$unsigned
    pval_FL_onesided[l] = 2*min(tmp1partial$pos,tmp1partial$neg)
  }
  
  ## PRegs t-test
  pval_PRegs_twosided_tstat = rep(NA, ncol(dat$y))
  pval_PRegs_onesided_tstat = rep(NA, ncol(dat$y))
  for(l in 1:ncol(dat$y)){
    tmp1partial <- PREG(y=dat$y[,l], x=dat$x, z=dat$Z, add_intercept = T, B = B, statType = "partial")$res
    pval_PRegs_twosided_tstat[l] = tmp1partial$unsigned
    pval_PRegs_onesided_tstat[l] = 2*min(tmp1partial$pos,tmp1partial$neg)
  }
  
  ## PRegs coef
  pval_PRegs_twosided_coef = rep(NA, ncol(dat$y))
  pval_PRegs_onesided_coef = rep(NA, ncol(dat$y))
  for(l in 1:ncol(dat$y)){
    tmp1coef <- PREG(y=dat$y[,l], x=dat$x, z=dat$Z, add_intercept = T, B = B, statType = "coef")$res
    pval_PRegs_twosided_coef[l] =  tmp1coef$unsigned
    pval_PRegs_onesided_coef[l] = 2*min(tmp1coef$pos,tmp1coef$neg)
  }
  
  ##cyclic permutations
  ### genetic alg ordering
  rounds1 <- 100
  rounds2 <- 900
  M = NULL
  X0 = cbind(dat$x, dat$Z[,-ncol(dat$Z)])
  X = cbind(dat$x, dat$Z)
  m = 20
  n = nrow(X)
  res1 <- find_eta_GA(X0, m, testinds = 1,
                      popSize = 10, rounds = rounds1,
                      M = M)
  res2 <- find_eta_GA(X0, m, ga_obj = res1$ga_obj,
                      testinds = 1,
                      popSize = 10, rounds = rounds2,
                      M = M)
  orderings = cbind(cbind(sample(1:n), res1$ordering),res2$ordering)
  
  pval_CPT_mat = matrix(NA, ncol  = 3, nrow = 3)
  alphas = c(0.05, 0.01, 0.001)
  l=1
  for(k in 1:length(alphas)){
    for(j in 1:ncol(orderings)){
      out0 = CPT(dat$y[,l], X, ordering = orderings[,j],testinds = 1,
                 alpha = alphas[k],
                 returnCI = FALSE)
      pval_CPT_mat[k,j] = out0$pval[1]
      if(out0$O<=0){
        pval_CPT_mat[k,j] = NA
      }
    }
  }
  pvals_CPT_array = array(NA, dim = c(ncol(dat$y), 3, 3))
  for(k in 1:length(alphas)){
    for(j in 1:ncol(orderings)){
      if(!is.na(pval_CPT_mat[k,j])){
        for(l in 1:ncol(dat$y)){
          out0 = CPT(dat$y[,l], X, ordering = orderings[,j],testinds = 1,
                     alpha = alphas[k],
                     returnCI = FALSE)
          pvals_CPT_array[l,k,j] =  out0$pval[1]
          if(out0$O==0){
            pvals_CPT_array[k,j] = NA
          }
        }
      }
    }
  }
  ##simple permutations
  pval_perm_twosided = rep(NA, ncol(dat$y))
  pval_perm_onesided = rep(NA, ncol(dat$y))
  for(l in 1:ncol(dat$y)){
    tmp1partial <- naivePermutation(y=dat$y[,l], x=dat$x, z=dat$Z, add_intercept = T, B = B)$res
    pval_perm_twosided[l] =  tmp1partial$unsigned
    pval_perm_onesided[l] = 2*min(tmp1partial$pos,tmp1partial$neg)
  }
  pvals_list = list()
  
  pvals_list[["twosided"]] = list()
  pvals_list[["twosided"]] = data.frame("ttest" = pval_ttest, "FL" = pval_FL_twosided,
                                        "PRegs" = pval_PRegs_twosided_tstat,
                                        "PRegsCoef" = pval_PRegs_twosided_coef,
                                        "perm" = pval_perm_twosided,
                                        "cyclicStrong05" = pvals_CPT_array[,1,3],
                                        "cyclicWeak05" = pvals_CPT_array[,1,2],
                                        "cyclicRandom05" = pvals_CPT_array[,1,1],
                                        "cyclicStrong01" = pvals_CPT_array[,2,3],
                                        "cyclicWeak01" = pvals_CPT_array[,2,2],
                                        "cyclicRandom01" = pvals_CPT_array[,2,1]
  )
  pvals_list[["onesided"]] = data.frame("ttest" = pval_ttest, "FL" = pval_FL_onesided,
                                        "PRegs" = pval_PRegs_onesided_tstat,
                                        "PRegsCoef" = pval_PRegs_onesided_coef,
                                        "perm" = pval_perm_onesided )
  apply( pvals_list[["twosided"]]<=0.05,2,mean)
  apply( pvals_list[["twosided"]]<=0.01,2,mean)
  apply( pvals_list[["twosided"]]<=0.001,2,mean)
  
  apply( pvals_list[["onesided"]]<=0.05,2,mean)
  apply( pvals_list[["onesided"]]<=0.01,2,mean)
  apply( pvals_list[["onesided"]]<=0.001,2,mean)
  return(pvals_list)
  
}
