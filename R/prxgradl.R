#' Penalized likelihood estimation of Cholesky factor
#' 
#' Solve the following optimization problem
#' \deqn{\hat{L} = \arg \min_{L} 2\log(det(L)) + tr(L^{-t}L^{-1} \Sigma)} + ||L||_1,off
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
#' @param job if 0 no additional zeros will be imposed
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
                        lambda = 0, job = 1){
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
cholpath <- function(Sigma, lambdas = NULL, L =  diag(nrow(Sigma)),
                    eps = 1e-8, maxIter = 1000){
  if (is.null(lambdas)) {
    lambdas = seq(0, max(diag(Sigma)), length = 10)
  }
  results <- list()
  for (i in 1:length(lambdas)){
    results[[i]] <- prxgradchol(Sigma, L, eps, 
                              maxIter = maxIter, lambda = lambdas[i])
    #L <- results[[i]]$L
  }
  return(results)
}