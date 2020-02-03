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
	p <- 10; d <- 0.5; N <- 10000
	
	B <- rlower(p = p, d = d)
	L <- diag(p) + B
 	S <- L %*% t(L)
 	data <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = S)
 	
 	reg <- fit_ldl_cov(amat = B, data = data)
 	L_chol <- cholfromldl(L = reg$L, D = reg$D)
 	L_lik <- prxgradchol(Sigma = cov(data), L = L_chol,
 										 lambda = 0, 
 										 maxIter = 1000, eps = 1e-15)
 	expect_equal(norm(L_lik$L - L_chol), 0, tolerance = 1e-2)
 	
 	ldl <- ldlfromchol(L = L_lik$L)
 	expect_equal(norm(reg$L - ldl$L), 0, tolerance = 1e-2)
 	expect_equal(norm(reg$D - ldl$D), 0, tolerance = 1e-2)
})
