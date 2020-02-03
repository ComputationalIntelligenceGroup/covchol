#' Returns LDL factorization from lower Cholesky factor
#' 
#' @param L Lower Cholesky factor
#' 
#' @return a list with the L and D factors (D being a vector)
#' @export
ldlfromchol <- function(L) {
	L_new <- L %*% diag(1/diag(L))
	D <- diag(L)^2
	
	return(list(L = L_new, D = D))
}

#' Returns Cholesky factorization from ldl
#' 
#' @param L L factor in LDL^t decomposition
#' @param D D factor (diagonal vector) in LDL^t decomposition
#' 
#' @return The lower Cholesky factor
#' @export
cholfromldl <- function(L,  D) {
	L_new <- L %*% diag(sqrt(D))
	
	return(L_new)
}

#' Simulates a lower triangular matrix, possibly with a sparsity pattern
#' 
#' @param p Dimension of the matrix to be simulated
#' @param d Density degree
#' 
#' @details Each entry in the lower triangle of the matrix is first assigned
#' the result of a Bernouilli with parameter `d`. If a `1` was obtained, then
#' the entry is assigned the result of a Uniform distribution on `[0.1, 1]`.
#' 
#' @return The simulated lower triangular matrix
#' @export
rlower <- function(p, d) {
	B <- matrix(nrow = p, ncol = p, data = 0)
	lower_tri_size <- p*(p - 1)/2
	B[lower.tri(B)] <- stats::rbinom(n = lower_tri_size, size = 1, prob = d) * 
		stats::runif(n = lower_tri_size, min = 0.1, max = 1)
	
	return(B)
}

#' Calculates the inverse of the Cholesky factor via the
#' regression coefficients formula
#' 
#' @param coeffs Lower-triangular matrix of regression coefficients
#' 
#' @return Matrix inverse of the Cholesky factor
#' @export
chol_inv <- function(coeffs) {
	o <- diag(nrow(coeffs))
	dimnames(o) <- list(rownames(coeffs), colnames(coeffs))
	p <- ncol(o)

	for (i in 2:p) {
		o[i, i - 1] <- coeffs[i, i - 1]
	}

	for (j in 1:(p - 2)) {
		for (i in (j + 2):p) {
			o[i, j] <- coeffs[i, j] +
				coeffs[i, (j + 1):(i - 1)] %*%
					o[(j + 1):(i - 1), j]
		}
	}

	return(o)
}

#' Fit the ldl factor of the inverse covariance / concentration matrix 
#' with linear regression
#'
#' @param amat Adjacency matrix marking the sparsity of the factor
#' @param data Data for fitting
#'
#' @return A list with the ldl factors
#' @export
fit_ldl_conc <- function(amat, data) {

	D <- numeric(nrow(amat))
	B <- amat
	for (i in seq(nrow(B), from = 2)) {
		parents <- which(B[i, ] != 0)
		if (length(parents) > 0) {
			model <- summary(stats::lm(data[, i] ~ data[, parents]))
			B[i, parents] <- - model$coefficients[- 1, 1]
			D[i] <- model$sigma^2
		} else {
			D[i] <- stats::var(data[, i])
		}
	}
	D[1] <- stats::var(data[, 1])
	diag(B) <- 1
	
	return(list(L = B, D = D))
}

#' Fit the ldl factor of the covariance matrix with linear regression
#' 
#' @param amat Adjacency matrix marking the sparsity of the factor
#' @param data Data for fitting
#' 
#' @return A list with the ldl factors
#' @export
fit_ldl_cov <- function(amat, data) {

	D <- numeric(ncol(amat))
	B <- amat
	for (i in seq(nrow(B), from = 2)) {
		parents <- which(B[i, ] != 0)
		for (j in parents) {
			model <- stats::lm(data[, i] ~ data[, 1:j])
			B[i, j] <- model$coefficients[j + 1]
		}
		model <- summary(stats::lm(data[, i] ~ data[, seq(i - 1)]))
		D[i] <- model$sigma^2
	}
	D[1] <- stats::var(data[, 1])
	diag(B) <- 1
	
	return(list(L = B, D = D))
}

