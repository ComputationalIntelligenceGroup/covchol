#' Minus Log-Likelihood for Gaussian model
#'
#' @param P the inverse of the covariance matrix
#' @param S the empirical covariance
#'
#' @importFrom stats var
#' @export
mll <- function(P, S){
  -determinant(P, logarithm = TRUE)$modulus + sum(S * P)
}

#' Proximal gradient sparse Cholesky factor
#'
#'
#' @param Sigma the empirical covariance matrix
#' @param L an initial guess for L
#' @param D the noise matrix 
#' @param eps stopping criteria
#' @param alpha parameter for line search
#' @param beta parameter for line search
#' @param maxIter maximum number of iterations
#' @param trace if >0 print info (1:termination message only, 
#'               2:messages at each step)
#' @param lambda penalization coefficient
#' @param r logical, if TRUE the set of non-zero parameters will be updated 
#'          at every iteration
#' @param h logical if TRUE only the non-zero entries of L will be updated
#'
#' @return the estimated L matrix 
#' @export
proxgradL <- function(Sigma, L, D = diag(ncol(Sigma)), eps =  1e-4,
                      alpha = 0.2, beta = 0.5,
                      maxIter = 100, trace = 0,
                      lambda = 0, r = FALSE, h = TRUE){
  p <- ncol(Sigma)
  ixd <- 0:(p - 1) * (p) + 1:p ##index of diagonal elements
  if (h) ix <- (1: p^2 )[L != 0]
  else ix <- 1:(p^2)
  ixnd <- ix[! ix %in% ixd ] ## not diagonal elements
  a <- Inf
  n <- 0
  S <- L %*% D %*% t(L)
  P <- solve(S)
  tk <- 1
  u <- rep(0, length(ix))
  f <- mll(P = P, S = Sigma)
  while (a > eps && n < maxIter) {
    if (r){
      ix <- (1: p^2 )[L != 0]
      ixnd <- ix[! ix %in% ixd ] ## not diagonal elements
    }
    n <- n + 1
    gradSigma <- P %*% Sigma %*% P - P 

    u <- 2 * gradSigma %*% D %*% t(L)
    diag(u) <- 0
    Lold <- L
    
    #### Beck and Teboulle line search
    f <- mll(P, Sigma) + lambda * sum(abs(Lold[ixnd]))
    fnew <- Inf
    
    alph <- alpha
    while ( fnew   > f - sum(u[ixnd] * (L[ixnd] - Lold[ixnd]) ) +
            sum((L[ixnd] - Lold[ixnd]) ^ 2) / (2* alph) || fnew > f) {
      
      L[ixnd] <- Lold[ixnd] + alph * u[ixnd]
      
      ### soft thres
      L[ixnd] <- sign(L[ixnd]) * (abs(L[ixnd]) - alph * lambda)
      L[ixnd][abs(L[ixnd]) < (alph * lambda)] <- 0
      
      ### new S
      S <- L %*% D %*% t(L)
      P <- solve(S)
      fnew <- mll(P, Sigma) + lambda * sum(abs(L[ixnd]))
      alph <- alph * beta
    }
    
    a <- (f - fnew ) / (abs(f))      
    if (trace > 1){
      message("Iteration: ", n, " ||diff||:", signif(a),
              " alpha:", alph / beta, "||L||_0:", sum(L!=0))
    }
  }
  if (trace > 0){
    message("Stop after ", n, " iterations, with ||diff||=", signif(a))
  }
  return(L)
}



#' Path of L estimates (log-lik)
#' 
#' @param Sigma empirical covariance matrix
#' @param lambdas increasing sequence of lambdas
#' @param D the known noise matrix 
#' @param L0 initial L matrix
#' @param eps convergence criteria
#' @param maxIter maximum iterations for each proximal gradient
#' @param r as in proxgradL
#' @param h as in proxgradL
#' @export
llpathL <- function(Sigma, lambdas = NULL, 
                    D = diag(nrow(Sigma)),
                    L0 = NULL,
                    eps = 1e-8, maxIter = 1000, 
                    r = FALSE, h = TRUE){
  if (is.null(lambdas)) {
    lambdas = seq(0, max(diag(Sigma)), length = 10)
  }
  results <- list()
  if (is.null(L0)){
    L0 <- t(chol(Sigma))
  }
  for (i in 1:length(lambdas)){
    results[[i]] <- proxgradL(Sigma, L0, D, eps, beta = 0.5, 
                                maxIter = maxIter, lambda = lambdas[i], r = r, 
                              h = h )
    L0 <- results[[i]]
    
  }
  return(results)
} 
