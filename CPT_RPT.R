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
