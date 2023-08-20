
##Permutation test
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

PREGS_CI_table = function(Cmat){
  tmat_collection = PREGS_CI_table_construct(Cmat)
  for(l in 1:length(tmat_collection)){
    tmat_collection[[l]] = data.frame(tmat_collection[[l]])
    colnames(tmat_collection[[l]]) = c("tvals", "scount", "ucount", "gamma", "gamma_plus",
                                       "pval", "pval_plus", "pval_low", "pval_high")
  }
  return(tmat_collection)
}
# ## PRegs CI construction
# PREGS_CI_table = function(Cmat){
#   ###
#   C1mat = Cmat[,1,]
#   C2mat = Cmat[,2,]
#   C3mat = Cmat[,3,]
#   C4mat = Cmat[,4,]
#   s = (C2mat - sqrt(C2mat^2 - C1mat * (C3mat - C4mat)))/C1mat
#   u = (C2mat + sqrt(C2mat^2 - C1mat * (C3mat - C4mat)))/C1mat
#   M = dim(Cmat)[3]
#   B = dim(Cmat)[1]
#   tmat_collection = list()
#   for(l in 1:M){
#     a1 = which(C1mat[,l] == 0 &(C3mat[,l] == C4mat[,l]))
#     a2 = which(C1mat[,l] == 0 &(C3mat[,l] < C4mat[,l]) )
#     a = setdiff(1:B, c(a1, a2))
#     sl = s[a,l]
#     ul = u[a,l]
#     sl = sort(sl)
#     ul = sort(ul)
#     tvals = unique(sort(c(sl, ul)))
#     tab_su = data.frame(tvals = tvals, scount = 0, ucount=0, gamma = 0, gamma_plus=0)
#     j = 1
#     for(i in 1:length(sl)){
#       while(sl[i] > tvals[j]){
#         j = j+1
#       }
#       tab_su$scount[j] = tab_su$scount[j]+1
#     }
#     j = 1
#     for(i in 1:length(ul)){
#       while(ul[i] > tvals[j]){
#         j = j+1
#       }
#       tab_su$ucount[j] = tab_su$ucount[j]+1
#     }
#     tab_su$gamma[1] = 0.5 * (tab_su$scount[1] - tab_su$ucount[1])
#     for(i in 2:nrow(tab_su)){
#       tab_su$gamma_plus[i-1] = tab_su$gamma[i-1]+0.5 * (tab_su$scount[i-1] - tab_su$ucount[i-1])
#       tab_su$gamma[i] = tab_su$gamma_plus[i-1]+0.5 * (tab_su$scount[i] - tab_su$ucount[i])
#     }
#     tab_su$pval =  (tab_su$gamma+0.5*length(a1)+length(a2)+1)/(B+1)
#     tab_su$pval_plus =  (tab_su$gamma_plus+0.5*length(a1)+length(a2)+1)/(B+1)
#     tab_su$pval_low = ifelse(tab_su$pval_plus>tab_su$pval,tab_su$pval_plus,tab_su$pval)
#     tmp1 = c(0, tab_su$pval_plus[-length(tab_su$pval_plus)])
#     tab_su$pval_high = ifelse(tmp1>tab_su$pval,tmp1, tab_su$pval)
#     tmat_collection[[l]] = tab_su
#   }
#   return(tmat_collection)
# }

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



