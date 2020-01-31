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
	S <- solve(t(L) %*% L)
	data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)
	
	U_reg <- fit_chol_conc(amat = B, data = data)
	colnames(data) <- colnames(B) <- rownames(B) <- seq(p)
	U_ggm <- ggm::fitDag(amat = t(B), S = cov(data), n = N)
	
	expect_equal(norm(U_reg - U_ggm$Ahat), 0)
})

test_that("estimation of covariance is correct", {
	p <- 10; d <- 0.5; N <- 100
	
	L <- rlower(p = p, d = d)
	diag(L) <- 1
	S <- L %*% t(L)
	data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)
	
	L_reg <- fit_chol_cov(amat = L, data = data)
	L_lik <- proxgradL(Sigma = cov(data), L = L_reg, D = diag(p),
										 lambda = 0, h = TRUE, trace = 2, alpha = 0.02, beta = 0.05, 
										 maxIter = 100, eps = 1e-10)
	
	expect_equal(norm(L_reg - L_lik$L), 0)
})
