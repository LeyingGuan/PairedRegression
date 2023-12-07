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



