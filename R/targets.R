#' Calculate the Target function for Regression
#'
#' @param y A vector of observations
#' @param X A design matrix
#' @param beta A vector of regression parameters
#'
#' @return The mean square error (MSE) from the linear model
#' @export
#'
#'

target_fun <- function(y, X, beta) {
    return(sum((y - X %*% beta)^2) / length(y))
}

#' Calculate the Gradient for Regression
#'
#' @inheritParams gradient_fun_penalty
#' @inheritParams target_fun
#'
#' @return The gradient vector
#' @export
#'

gradient_fun <- function(y, tXy, tXX, beta) {
    return((tXX %*% beta - tXy) / length(y))
}


#' Calculate the Gradient for Penalized Regression
#'
#' @inheritParams target_fun
#' @param tXy A vector of the transpose of design matrix X times the output data y
#' @param tXX A matrix of the  transpose of design matrix X times the design matrix X
#' @param Q A penalty (precision) matrix
#'
#' @return The gradient vector
#' @export
#'

gradient_fun_penalty <- function(y, tXy, tXX, beta, Q) {
    return((tXX %*% beta - tXy + Q %*% beta) / length(y))
}
