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

