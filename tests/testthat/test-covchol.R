test_that("estimation of covariance is correct", {
  p <- 5; d <- 1; N <- 10000;
  
  L <- gmat::mh_u(N = 1, p = p, dag = gmat::rgraph(p = p, d = d, dag = TRUE))[, , 1][p:1, p:1]
  S <- L %*% t(L)
  X <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)
  
  out <- .Fortran("PRXGRD",as.integer(ncol(X)), as.double(S), 
                  as.double(L), as.double(0), as.double(1e-15),
                  as.double(0.5), as.integer(1000),
                  PACKAGE = "covchol")
  names(out) <- c("N","Sigma", "L", "lambda", "diff", 
                  "objective", "iter")
  out$L <- matrix(nrow = out$N, out$L)
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  
  expect_equal(norm(out$Sigma - S, type = "F"), 0)
})


test_that("estimation of covariance with noise is correct", {
  p <- 5; d <- 1; N <- 10000;
  
  L <- gmat::mh_u(N = 1, p = p, dag = gmat::rgraph(p = p, d = d, dag = TRUE))[, , 1][p:1, p:1]
  S <- L %*% t(L)
  X <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)
  
  out <- .Fortran("PRXGRD",as.integer(ncol(X)), as.double(S), 
                  as.double(L + diag(runif(p, 0.01, 0.1))), as.double(0), as.double(1e-15),
                  as.double(0.5), as.integer(1000),
                  PACKAGE = "covchol")
  names(out) <- c("N","Sigma", "L", "lambda", "diff", 
                  "objective", "iter")
  out$L <- matrix(nrow = out$N, out$L)
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  
  expect_equal(norm(out$Sigma - S, type = "F"), 0)
})

test_that("optimization problem is correct lambda = 0 start from optimum", {
	p <- 5; d <- 1; 
	
 	S <- gmat::chol_mh(1, p = p,d = 1)[,,1]
  L <- t(chol(S))
 	out <- .Fortran("PRXGRD",as.integer(p), as.double(S), 
 									as.double(L), as.double(0), as.double(1e-6),
 									as.double(0.5), as.integer(100),
 									PACKAGE = "covchol")
 	names(out) <- c("N","Sigma", "L", "lambda", "diff", 
 									"objective", "iter")
 	out$L <- matrix(nrow = out$N, out$L)
 	out$Sigma <- matrix(nrow = out$N, out$Sigma)
 	
 	expect_equal(norm(out$Sigma - S, type = "F"), 0)
})


test_that("optimization problem is correct lambda = 0 start from optimum + noise", {
  p <- 5; d <- 1; 
  
  S <- gmat::chol_mh(1, p = p,d = 1)[,,1]
  L <- t(chol(S))
  L[lower.tri(L)] <- L[lower.tri(L)] + runif(sum(lower.tri(L))) 
  out <- .Fortran("PRXGRD",as.integer(p), as.double(S), 
                  as.double(L), as.double(0), as.double(1e-6),
                  as.double(0.5), as.integer(100),
                  PACKAGE = "covchol")
  names(out) <- c("N","Sigma", "L", "lambda", "diff", 
                  "objective", "iter")
  out$L <- matrix(nrow = out$N, out$L)
  out$Sigma <- matrix(nrow = out$N, out$Sigma)
  
  expect_equal(norm(out$Sigma - S, type = "F"), 0)
})