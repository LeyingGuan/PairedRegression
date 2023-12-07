library(plyr)
##Permutation test'
#' @useDynLib PREGS
#' @importFrom Rcpp sourceCpp
#' @export
PERMtest = function(y, x, z, add_intercept = T, B = 1E4){
  n = nrow(z)
  z = cbind(z, rep(1, n))
  if(is.null(dim(y))){
    y = matrix(y, ncol = 1)
  }
  out = permutation_vanilla(x=x, Y = y, Z=z, B=B)
  return(list(res = out))
}

## FL test
#' @useDynLib PREGS
#' @importFrom Rcpp sourceCpp
#' @export
FLtest = function(y, x, z, add_intercept = T, B = 1E4, statType = "coef"){
  n = nrow(z)
  z = cbind(z, rep(1, n))
  if(is.null(dim(y))){
    y = matrix(y, ncol = 1)
  }
  out = permutation_FL(x=x, Y = y, Z=z, B=B)
  return(list(res = out))
  
}
## RPT
generatePermutations <- function(n, K){
  ind <- head(sample(n), n - n %% (K+1))
  mx <- matrix(ind, nrow=K+1)
  perm_ind <- matrix(0, K, n)
  for (r in 1:K){
    perm_ind[r, ] <- plyr::mapvalues(1:n, mx, mx[c((r+1):(K+1), 1:r), ])
  }
  return(perm_ind)
}


## RPT test
#' @export
RPT <- function(X, Y, Z, K=99){
  n <- dim(X)[1]; p <- dim(X)[2]
  if(!is.null(dim(Y))){
    stat1 <- stat2 <- matrix(0,ncol=K, nrow = ncol(Y))
      rep(0, K)
  }else{
    stat1 <- stat2 <- rep(0, K)
  }
  
  #stat1[1] <- stat2[1] <- sum(ehat * epshat)
  perm_ind <- generatePermutations(n, K)
  
  for (r in 1:K){
    idx <- perm_ind[r, ]  # P is a perm matrix with P[i, idx[i]] = 1
    #V <- U[idx, ]  # V = PU
    
    #tmp <- cbind(X, X[idx,])
    #VtildeVtilde <- diag(n) - tmp %*% solve(t(tmp) %*% tmp, t(tmp))
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
    #printPercentage(r, K)
  }
  if(is.null(dim(stat1))){
    stat1 = matrix(stat1,nrow = 1)
    stat2 = matrix(stat2,nrow = 1)
  }
  pval <-sapply(1:nrow(stat1), function(i){ (sum(abs(stat2[i,]) >= abs(stat1[i,])) + 1) / (K+1)})
  pval_theory <-sapply(1:nrow(stat1), function(i){ (sum(abs(stat2[i,]) >= min(abs(stat1[i,]))) + 1) / (K+1)})
  return(list(pval=pval, pval_theory=pval_theory))
}

## PRegs
#' @useDynLib PREGS
#' @importFrom Rcpp sourceCpp
#' @export
PREGtest = function(y, x, z, add_intercept = T, B = 1E4){
  n = nrow(z)
  z = cbind(z, rep(1, n))
  if(is.null(dim(y))){
    y = matrix(y, ncol = 1)
  }
  out = permutation_PREGSjoint(x=x, Y = y, Z=z, B=B)
  return(list(res = out))
}

## PRegs CI table
#' @useDynLib PREGS
#' @importFrom Rcpp sourceCpp
#' @export
PREGS_CI_table = function(Cmat){
  tmat_collection = PREGS_CI_table_construct(Cmat)
  for(l in 1:length(tmat_collection)){
    tmat_collection[[l]] = data.frame(tmat_collection[[l]])
    colnames(tmat_collection[[l]]) = c("tvals", "scount", "ucount", "gamma", "gamma_plus",
                                       "pval", "pval_plus", "pval_low", "pval_high")
  }
  return(tmat_collection)
}

## PREG CI contruction
#' @export
PREG_CI_construct= function(tmat_collection, alpha = 0.05){
  CI = data.frame(matrix(NA, ncol = 2, nrow = length(tmat_collection)))
  colnames(CI) = c("beta_min", "beta_max")
  for(l in 1:length(tmat_collection)){
    tab_su = tmat_collection[[l]]
    ll = which(tab_su$pval_low> alpha)
    if(length(ll) > 0){
      id_min = min(which(tab_su$pval_low> alpha))
      id_max = max(which(tab_su$pval_high> alpha))
      CI[l,1] = tab_su$tvals[id_min]
      CI[l,2] = tab_su$tvals[id_max]
    }
  }
  return(CI)
}
# tmat_collection = PREGS_CI_table(Cmat)
# CI = PREGci_cut(tmat_collection, alpha = 0.05)
# 1 - mean(CI[,2]> 0 & CI[,1]<0)



