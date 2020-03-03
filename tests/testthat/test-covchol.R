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