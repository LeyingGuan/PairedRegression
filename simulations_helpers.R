library(PREGS)
library(MASS)

data_gen = function(n=100, p = 2, M = 2000, design = "AnovaBalance", noise = "gaussian",seed = 1, m =5){
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
sim_comparisons_singleSetting = function(dat, run_CPT=TRUE, B = 2000){
  
  pvals_collection = data.frame(matrix(NA, ncol = 5+1, nrow = ncol(dat$y)))
  colnames(pvals_collection) = c("Ftest", "PERMtest", 
                                 "FLtest", "PREGStest","PREGStestII",
                                 "CPT")
  ###t test
  fitted_ttest = summary(lm(dat$y~dat$x+dat$Z))
  pvals_collection[,1] = as.vector(sapply(fitted_ttest, function(z) z$coefficients[2,4]))
  
  ##simple permutations
  pvals_collection[,2] = PERMtest(x= dat$x, y = dat$y, z=dat$Z, B = B)$res$unsigned
  ###FL-test (F)
  pvals_collection[,3] =  FLtest(x= dat$x, y = dat$y, z=dat$Z, B = B)$res$unsigned
  
  ## PRegs
  pvals_collection[,4] = PREGtest(x= dat$x, y = dat$y, z=dat$Z, B = B, mode = "separate")$res$unsigned
  pvals_collection[,5] = PREGtest(x= dat$x, y = dat$y, z=dat$Z, B = B, mode = "joint")$res$unsigned


  ##cyclic permutations
  ### genetic alg ordering
  if( run_CPT){
    rounds <- 100 + 900; m = 19;M = NULL
    X0 = cbind(dat$x, dat$Z[,-ncol(dat$Z)])
    X = cbind(dat$x, dat$Z)
    res1 <- find_eta_GA(X0, m, testinds = 1, popSize = 10, rounds = rounds,M = M)
    ordering = res1$ordering
    for(l in 1:ncol(dat$y)){
      out0 = CPT(dat$y[,l], X, ordering = ordering,testinds = 1,
                 alpha = 0.05,
                 returnCI = FALSE)
      pvals_collection[l,6]=  out0$pval[1]
    }
  }
  # apply(pvals_collection<=0.05,2,mean)
  # apply(pvals_collection<=0.01,2,mean)
  # apply(pvals_collection<=0.001,2,mean)

  return(pvals_collection)
  
}
