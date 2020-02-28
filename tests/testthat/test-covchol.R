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
	p <- 10; d <- 1; N <- 10000;
	
	L <- rlower(p = p, d = d)
	diag(L) <- runif(p, 0.1, 1)
 	S <- L %*% t(L)
 	X <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)

 	L_scale <- diag(1/sqrt(diag(cov(X)))) %*% L
	res <- prxgradchol(X, L = L_scale, lambda = 0)
	
 	expect_equal(norm(res$Sigma - S, type = "F"), 0, tolerance = 1e-2)
})
