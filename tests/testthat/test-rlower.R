test_that("random lower triangular matrix is lower triangular", {
  p <- 10; d <- 0.5
  
  B <- rlower(p = p, d = d)
  
  expect_equal(sum(abs(B[upper.tri(B, diag = TRUE)])), 0)
})
