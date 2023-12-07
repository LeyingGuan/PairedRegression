library("MASS")
library("RSpectra")
## devtools::install_url("https://cran.r-project.org/src/contrib/Archive/gaoptim/gaoptim_1.1.tar.gz")
library("gaoptim")

## Calculate the 2-norm of a vector
norm2 <- function(x){
  sqrt(sum(x^2))
}

## Shift a vector 1:n to the left.
left_shift <- function(n, m){
  quotient <- floor(n / (m + 1))
  tempn <- (m + 1) * quotient
  remainder <- n - tempn
  sapply(0:m, function(i){
    nshift <- i * quotient
    shiftseq <- 1:tempn + nshift - c(rep(0, tempn - nshift), rep(tempn, nshift))
    if (remainder > 0){
      shiftseq <- c(shiftseq, tempn + 1:remainder)
    }
    shiftseq
  })
}

## Transform the design matrix so that the linear hypothesis
## R\beta = 0 becomes equivalent to \beta_1 = ... = \beta_r
## = 0
#' @export
preprocess_X <- function(X, R = NULL, testinds = NULL){
  if (!is.null(R)){
    r <- nrow(R)
    p <- ncol(R)
    obj <- svd(t(R))
    neigs <- sum(abs(obj$d) > 1e-10)
    U <- obj$u[, 1:neigs]
    V <- svd(diag(r) - U %*% t(U))$u[, 1:(p - neigs)]
    X <- X %*% cbind(U, V)
    ntestinds <- neigs
  } else if (!is.null(testinds)){
    testinds <- sort(testinds)
    ntestinds <- length(testinds)
    if (!all(testinds == 1:ntestinds)){
      X <- cbind(X[, testinds], X[, -testinds])
    }
  } else {
    stop("Either R or testinds should be given.")
  }
  return(list(X = X, ntestinds = ntestinds))
}

## Find optimal eta given X
#' @export
solve_optim_eta <- function(X, m, ntestinds,
                            M = NULL){
  n <- nrow(X)
  Piinds <- left_shift(n, m)
  BX <- lapply(1:m, function(i){
    X[Piinds[, i], ] - X[Piinds[, m + 1], ]})
  BX <- do.call(cbind, BX)
  mod <- lm(BX[, 1:ntestinds] ~ BX[, -(1:ntestinds)] + 0)
  reg_resid <- resid(mod)
  if (ntestinds == 1){
    Ostar <- norm2(reg_resid)
    etastar <- reg_resid / Ostar
  } else {
    ## if (is.null(M)){
    ##     M <- t(X[, 1:ntestinds]) %*% X[, 1:ntestinds]
    ## }
    MrX <- reg_resid %*% M %*% t(reg_resid)
    obj <- RSpectra::eigs_sym(MrX, 1)
    Ostar <- obj$values
    etastar <- obj$vectors / norm2(obj$vectors)
  }
  return(list(etastar = etastar, Ostar = Ostar))
}

## Find optimal ordering using stochastic search
#' @export
find_eta_SS <- function(X, m,
                        R = NULL, testinds = NULL,
                        M = NULL,
                        ntry = 1000){
  n <- nrow(X)    
  best_Ostar <- 0
  best_ordering <- 1:n
  temp <- preprocess_X(X, R, testinds)
  X <- temp$X
  ntestinds <- temp$ntestinds
  for (i in 1:ntry){
    ordering <- sample(n, n)
    Ostar <- solve_optim_eta(X[ordering, ], m, ntestinds, M)$Ostar
    if (Ostar > best_Ostar){
      best_ordering <- ordering
      best_Ostar <- Ostar
    }
  }
  return(best_ordering)
}

## Find optimal ordering using genetic algorithm
#' @export
find_eta_GA <- function(X, m,
                        ga_obj = NULL,
                        R = NULL, testinds = NULL,
                        M = NULL,
                        rounds = 100,
                        popSize = 10, ...){
  if (is.null(ga_obj)){
    n <- nrow(X)
    temp <- preprocess_X(X, R, testinds)
    X <- temp$X
    ntestinds <- temp$ntestinds
    objective_fun <- function(ordering){
      solve_optim_eta(X[ordering, ], m, ntestinds, M)$Ostar
    }
    ga_obj <- gaoptim::GAPerm(objective_fun, n, popSize = popSize)
  }
  ga_obj$evolve(rounds)
  best_ordering <- ga_obj$bestIndividual()
  return(list(ordering = best_ordering,
              ga_obj = ga_obj))
}



## Cyclic Permutation Test
#' @export
CPT <- function(y, X,
                ordering = NULL,
                R = NULL, testinds = NULL,
                M = NULL,
                alpha = 0.05,
                m = ceiling(1 / alpha) - 1,
                returnCI = FALSE,
                ...){
  n <- nrow(X)
  if (is.null(ordering)){
    ordering <- 1:n
    warning("ordering is missing. Use 1:nrow(X) as default. May try to find a better ordering by find_eta_GA")
  }
  temp <- preprocess_X(X, R, testinds)
  X <- temp$X
  ntestinds <- temp$ntestinds
  
  y <- y[ordering]
  X <- X[ordering, ]
  tmp = solve_optim_eta(X, m, ntestinds)
  eta <- tmp$etastar
  Piinds <- left_shift(n, m)
  stats <- apply(Piinds, 2, function(inds){
    sum(y[inds] * eta)
  })
  stats <- abs(stats - median(stats))
  
  ## Two-sided p-value
  pval <- rank(-stats, ties.method = "max") / (m + 1)
  return(list(pval = pval,
              eta = eta,
              ordering = ordering,
              O =tmp$Ostar))
}
