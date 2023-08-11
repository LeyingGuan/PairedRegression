
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
PREGtest = function(y, x, z, add_intercept = T, B = 1E4, mode = "separate"){
  n = nrow(z)
  z = cbind(z, rep(1, n))
  if(is.null(dim(y))){
    y = matrix(y, ncol = 1)
  }
  if(mode == "separate"){
    out = permutation_PREGSseparate(x=x, Y = y, Z=z, B=B)
  }else{
    out = permutation_PREGSjoint(x=x, Y = y, Z=z, B=B)
  }

  return(list(res = out))
}


