#' Simulate a B-spline mixed model
#'
#' @param N The number of observations
#' @param df The B-spline degrees of freedom
#' @param sigma The residual error standard deviation
#' @param sigma_alpha The spline random effect standard deviation
#'
#' @return A list containing the simulated data vector y, the simulated design matrix X,
#' the simulated vector of fixed effects beta, the simulated B-spline basis Z,
#' the simulated random effects alpha, and the residual errors epsilon
#' @export
#'
#' @import splines
#'
#'
sim_example_splines <- function(N = 500, df = 8, sigma = 0.5, sigma_alpha = 1.5) {
    X <- cbind(1, rnorm(N))
    Z <- bs(seq(0, 1, length.out = N), df = 8)
    beta <- rnorm(ncol(X))
    alpha <- rnorm(ncol(Z), 0, sigma_alpha)
    epsilon <- rnorm(N, 0, sigma)
    y <- X %*% beta + Z %*% alpha + epsilon
    return(list(y = y, X = X, beta = beta, Z = Z,
                alpha = alpha, epsilon = epsilon))
}
