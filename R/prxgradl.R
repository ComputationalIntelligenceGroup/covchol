#' Penalized likelihood estimation of Cholesky factor
#' 
#' Solve the following optimization problem
#' \deqn{\hat{L} = \arg \min_{L} 2\log(det(L)) + tr(L^{-1} \Sigma) L^{-t}} + ||L||_1,off
#' 
#' \code{cholpath} returns the path of regularized estimator on a sequence of 
#' \code{lambda} parameters
#' 
#' @param X data from which to obtain the path
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
prxgradchol <- function(X, L = diag(ncol(X)), eps =  1e-2,
                        alpha = 0.5, 
                        maxIter = 100, 
                        lambda = 0, job = 1) {
	Cor <- stats::cor(X)
  
	out <- .Fortran("PRXGRD",as.integer(ncol(X)), as.double(Cor), 
                  as.double(L), as.double(lambda), as.double(eps),
                  as.double(alpha), as.integer(maxIter),
                  PACKAGE = "covchol")
  names(out) <- c("N","Sigma", "L", "lambda", "diff", 
                  "objective", "iter")
  out$L <- matrix(nrow = out$N, out$L)

  # Return to covariance matrices
 	D_scale <- sqrt(diag(stats::cov(X)))
  out$L <- diag(D_scale) %*% out$L
  # compute Sigma
  out$Sigma <- tcrossprod(out$L)
  
  return(out)
}


#' @rdname prxgradchol
#' 
#' @param lambdas increasing sequence of lambdas
#' @export
cholpath <- function(X, lambdas = NULL, L =  diag(ncol(X)),
                    eps = 1e-8, maxIter = 1000){
  if (is.null(lambdas)) {
  	p <- ncol(X)
  	N <- nrow(X)
  	
  	if (p < N) {
  		lambda.min.exp <- -4
  	} else {
  		lambda.min.exp <- -2
  	}
  	
    lambdas <- 10^seq(0, lambda.min.exp, length = 100)
  }
  results <- list()
  for (i in 1:length(lambdas)){
    results[[i]] <- prxgradchol(X, L, eps, 
                              maxIter = maxIter, lambda = lambdas[i])
    L <- results[[i]]$L
  }
  return(results)
}