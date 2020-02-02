#' Penalized likelihood estimation of Cholesky factor
#' 
#' Optimize the  matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) using proximal gradient. 
#' \deqn{\hat{L} = \arg \min_{L} 2\log(det(L)) + tr(L^{-t}L^{-1} \Sigma)}
#' 
#' \code{cholpath} returns the path of regularized estimator on a sequence of 
#' \code{lambda} parameters
#' 
#' @param Sigma the empirical covariance matrix
#' @param L initial cholesky factor
#' @param eps convergence threshold for the proximal gradient
#' @param alpha line search rate
#' @param maxIter the maximum number of iterations
#' @param lambda penalization coefficient 
#' 
#' @return a list with the output of the optimization:
#' 
#' * \code{N}
#' * \code{L} the estimated L matrix
#' * \code{lambda} 
#' * \code{diff} the value of the last relative decrease
#' * \code{objective} the value of the objective function
#' * \code{iter} number of iterations
#' @useDynLib covchol
#' @export
prxgradchol <- function(Sigma, L, eps =  1e-2,
                        alpha = 0.5, 
                        maxIter = 100, 
                        lambda = 0){
  out <- .Fortran("PRXGRD",as.integer(ncol(Sigma)), as.double(Sigma), 
                  as.double(L), as.double(lambda), as.double(eps),
                  as.double(alpha), as.integer(maxIter),
                  PACKAGE = "covchol")
  names(out) <- c("N","Sigma", "L", "lambda", "diff", 
                  "objective", "iter")
  out$L <- matrix(nrow = out$N, out$L)
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  return(out)
}


#' @rdname prxgradchol
#' 
#' @param lambdas increasing sequence of lambdas
#' @export
cholpath <- function(Sigma, lambdas = NULL, 
                    eps = 1e-8, maxIter = 1000){
  L0 <- t(chol(Sigma))
  if (is.null(lambdas)) {
    lambdas = seq(0, max(diag(Sigma)), length = 10)
  }
  results <- list()
  for (i in 1:length(lambdas)){
    results[[i]] <- prxgradchol(Sigma, L0, eps, 
                              maxIter = maxIter, lambda = lambdas[i])
    L0 <- results[[i]]$L
  }
  return(results)
}