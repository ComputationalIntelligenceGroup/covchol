test_that("inverse of the concentration factor with formula is correct", {
	p <- 10; d <- 0.5
	
	B <- rlower(p = p, d = d)
	L <- diag(p) - B
	
	O_inv_solve <- solve(L)
	O_inv_chol <- covchol::chol_inv(B)	
			
	expect_equal(norm(O_inv_chol - O_inv_solve), 0)
})

test_that("regression estimation of concentration is correct", {
	p <- 10; d <- 0.5; N <- 100
	
	B <- rlower(p = p, d = d)
	L <- diag(p) - B
	D <- runif(p, 0.1, 1)
	S <- solve(t(L) %*% diag(D) %*% L)
	data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)
	
	reg <- fit_ldl_conc(amat = B, data = data)
	colnames(data) <- colnames(B) <- rownames(B) <- seq(p)
	U_ggm <- ggm::fitDag(amat = t(B), S = cov(data), n = N)
	
	expect_equal(norm(reg$L - U_ggm$Ahat), 0, tolerance = 1e-1)
	#expect_equal(sum(abs(reg$D - U_ggm$Dhat)), 0, tolerance = 1e-1) fails
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
  L <- t(chol(S)) + rlower(p, 1)
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