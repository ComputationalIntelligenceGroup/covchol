test_that("random lower triangular matrix is lower triangular", {
  p <- 10; d <- 0.5
  
  B <- rlower(p = p, d = d)
  
  expect_equal(sum(abs(B[upper.tri(B, diag = TRUE)])), 0)
})

test_that("ldl from Cholesky is correct", {
	p <- 5; d <- 0.5
	
	B <- rlower(p = p, d = d)
	L <- diag(p) + B
	D <- diag(runif(p, 0.1, 1))
	
	S <- L %*% D %*% t(L)
	L_chol <- t(chol(S))
	
	res <- ldlfromchol(L = L_chol)
	
	expect_equal(norm(res$L - L), 0)
	expect_equal(norm(res$D - D), 0)
})