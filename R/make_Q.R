#' Generate CAR precision matrix
#'
#' A function for setting up a conditional autoregressive (CAR) or simultaneous autoregressive (SAR) precision matrix for use as a prior in Bayesian models
#'
#' @param n_dims is a vector of length M that are the dimensions of the CAR/SAR matrix at each resolution
#' @param phi is a vector of length M with each element between -1 and 1 that defines the strength of the  autoregressive process. Typically this will be set to 1 for use as a prior in penalized Bayesian models
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE}(see Matrix package for details)
#' @param prec_model is a string that takes the values "CAR" or "SAR" and defines the graphical structure for the precision matrix.
#' @return a list of  \eqn{n \times n}{n x n} sparse spam matrices or Matrix matrices of class "dgCMatrix" (see Matrix package for details)
#' @importFrom igraph as_adjacency_matrix make_lattice
#' @importFrom Matrix Diagonal colSums
#' @importFrom spam spam
#'
#' @examples
#' n_dims <- c(4, 8)
#' phi <- c(0.8, 0.9)
#' Q <- make_Q(n_dims, phi)
#' ## plot the precision matrix structure at each resolution
#' layout(matrix(1:2, 1, 2))
#' spam::display(Q[[1]])
#' spam::display(Q[[2]])
#'
#' @export
make_Q <- function(n_dims, phi, use_spam = TRUE, prec_model = "CAR") {

    # if (any(phi < -1) || any(phi > 1) || any(is.na(phi)))
    #     stop("phi must be a numeric vector of length M with values between -1 and 1.")
    #
    # if (!is_integer(n_dims, length(n_dims)))
    #     stop("n_dims must be a vector of integers of length M.")
    #
    # if (length(n_dims) != length(phi))
    #     stop("n_dims and phi must both be vectors of length M.")
    #
    # if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam))
    #     stop("use_spam must be either TRUE or FALSE.")
    #
    # if (!(prec_model %in% c("CAR", "SAR")))
    #     stop('The only valid options for prec_model are "CAR" and "SAR".')


    M <- length(n_dims)
    Q <- vector(mode = "list", length = M)
    for (m in 1:M) {
        W <- as_adjacency_matrix(make_lattice(length = n_dims[[m]], dim = 2), sparse = TRUE)
        D <- Diagonal(x = colSums(W))
        if (prec_model == "CAR") {
            Q[[m]] <- D - phi[m] * W
        } else if (prec_model == "SAR") {
            B <- diag(nrow(W)) - phi * W %*% Diagonal(x = 1 / colSums(W))
            Q[[m]] <- t(B) %*% B
        }
        if (use_spam) {
            ## use the spam package for sparse matrices
            # Q[[m]] <- spam(c(as.matrix(Q[[m]])), nrow = n_dims[[m]]^2, ncol = n_dims[[m]]^2)
            Q[[m]] <- as.spam.dgCMatrix(Q[[m]])
        }
    }

    if (M == 1) {
        Q <- Q[[1]]
    }


    if (length(n_dims) > 1) {
        for (m in 1:M){
            Q[[m]] <- 2^(2*(m-1)) * Q[[m]]
        }
        if (use_spam) {
            Q <- do.call(bdiag.spam, Q)
        } else {
            Q <- do.call(bdiag, Q)
        }
    }

    return(Q)
}
