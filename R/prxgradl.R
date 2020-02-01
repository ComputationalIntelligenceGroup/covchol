#' Penalized likelihood estimation of Cholesky factor
#' 
#' Optimize the  matrix of a continuous Lyapunov 
#' Gaussian graphical model (CLGGM) using proximal gradient. 
#' \deqn{\hat{L} = \arg \min_L \log + \lambda||L||_1}
#' 
#' @param Sigma the observed covariance matrix
#' @param L an initial L matrix
#' @param eps convergence threshold for the proximal gradient
#' @param alpha Beck and Tabulle line search rate
#' @param maxIter the maximum number of iterations
#' @param lambda penalization coefficient 
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
prxgradl <- function(Sigma, L, eps =  1e-2,
                        alpha = 0.5, 
                        maxIter = 100, 
                        lambda = 0){
  out <- .Fortran("PRXGRDL",as.integer(ncol(Sigma)), as.double(Sigma), 
                  as.double(L), as.double(lambda), as.double(eps),
                  as.double(alpha), as.integer(maxIter),
                  PACKAGE = "covchol")
  names(out) <- c("N","Sigma", "L", "lambda", "diff", 
                  "objective", "iter")
  out$L <- matrix(nrow = out$N, out$L)
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  return(out)
}