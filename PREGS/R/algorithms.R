perm_prepare = function(y, x, z, add_intercept){
  q = ncol(z)
  p = ncol(x)
  n = ncol(z)
  if(add_intercept == T){
    Z1 = cbind(rep(1,n), z)
  }else{
    Z1 = z
  }
  Zsvd = svd(Z1)
<<<<<<< HEAD
  U = Zsvd$u[,Zsvd$d>=1E-10,drop = F] # Hz = U U^T
=======
  d = Zsvd$d
  U = Zsvd$u[,  d>=1E-10] # Hz = U U^T
>>>>>>> e55732742df7fa41a5b995396d7cf72b7cf6cee1
  yfitted = U%*%(t(U)%*%y)
  yresid = y - yfitted
  Xresid = x- U%*%(t(U)%*%x)
  return(list(U = U, yfitted = yfitted, yresid = yresid, Xresid = Xresid))
}

# Algorithms wrapper

## T test (lm regression model) not need to implement


##Permutation test
#' @export
naivePermutation = function(y, x, z, add_intercept = T, B = 1E4){
  n = nrow(z)
  out0 = perm_prepare(y, x, z, add_intercept=add_intercept)
  U = out0$U; yresid = out0$yresid; Xresid = out0$Xresid
  perm_idx = sapply(1:B, function(b) sample(1:n, n, replace = F))
  perm_idxC = perm_idx-1
  out =  permutation_simple_C(X= x, yresid=yresid, U=U, perm_idx=perm_idxC)
  return(list(res = out))
}
## FL test
#' @export
FL = function(y, x, z, add_intercept = T, B = 1E4, statType = "coef"){
  n = nrow(z)
  out0 = perm_prepare(y, x, z, add_intercept=add_intercept)
  U = out0$U; yfitted = out0$yfitted; yresid = out0$yresid; Xresid = out0$Xresid
  perm_idx = sapply(1:B, function(b) sample(1:n, n, replace = F))
  perm_idxC = perm_idx-1
  if(statType == "coef"){
    out =  permutation_FL_C(Xresid=Xresid, yfitted=yfitted,yresid=yresid, U=U, perm_idx=perm_idxC , type = "coef")
  }else{
    out =  permutation_FL_C(Xresid=Xresid, yfitted=yfitted,yresid=yresid, U=U, perm_idx=perm_idxC, type = "partial")
  }
  return(list(res = out))
  
}

## PRegs
#' @export
PREG = function(y, x, z, add_intercept = T, B = 1E4, statType = "coef"){
  n = nrow(z)
  out0 = perm_prepare(y, x, z, add_intercept=add_intercept)
  U = out0$U; yfitted = out0$yfitted; yresid = out0$yresid; Xresid = out0$Xresid
  perm_idx = sapply(1:B, function(b) sample(1:n, n, replace = F))
  perm_idxC = perm_idx-1
  if(statType == "coef"){
    out = permutation_conformal_C(Xresid=Xresid, yfitted=yfitted,yresid=yresid, U=U, perm_idx=perm_idxC, type = "coef")
  }else{
    out = permutation_conformal_C(Xresid=Xresid, yfitted=yfitted,yresid=yresid, U=U, perm_idx=perm_idxC, type = "partial")
  }
  return(list(res=out))
  
}

## CyclicPerm
