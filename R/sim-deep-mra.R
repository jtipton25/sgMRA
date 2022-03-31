#' Simulate deep MRA
#'
#' @param N The number of locations
#' @param M The number of resolutions in the MRA approximation to the Gaussian process
#' @param n_coarse_grid The number of grid locations in the coarsest resolution of the MRA approximation to the Gaussian process
#' @param n_layers The number of layers in the deep process
#' @param sigma The observation standard deviation
#' @param use_spam A logical value of whether to use the spam library (use_spam=TRUE) or the Matrix library (use_spam=FALSE) for sparse matrix calculations
#' @param ncores The number of cores to use for parallelization
#' @param nchunks The number of chunks to divide the distance calculation into. The default argument of NULL will use the same number of chunks as the number of cores.
#'
#' @return
#' @export
#'
#' @importFrom tidyr expand_grid
#'
sim_deep_mra <- function(N = 100^2, M = 1, n_coarse_grid = 20,
                         n_layers = 3, sigma = 0.25, use_spam = FALSE, ncores = 1L, nchunks = NULL) {

    # define the locations and grid
    locs <- expand_grid(x=seq(0, 1, length.out=sqrt(N)),
                        y=seq(0, 1, length.out=sqrt(N)))
    grid <- make_grid(locs, M = M, n_coarse_grid = n_coarse_grid)

    # initialize the MRA parameters
    MRA <- vector(mode = 'list', length = n_layers)
    Q <- vector(mode = 'list', length = n_layers)
    # CH <- vector(mode = 'list', length = n_layers)
    alpha_x <- vector(mode = 'list', length = n_layers-1)
    alpha_y <- vector(mode = 'list', length = n_layers-1)

    # initialize the bottom layer
    MRA[[n_layers]] <- eval_basis(locs, grid, use_spam = use_spam, ncores = ncores, nchunks = nchunks)
    Q[[n_layers]] <- make_Q(sqrt(MRA[[n_layers]]$n_dims),
                            phi=rep(0.9, length(MRA[[n_layers]]$n_dims)),
                            use_spam = use_spam)
    # if (!use_spam) {
    #     CH[[n_layers]] <- Cholesky(Q[[n_layers]])
    # }

    # if (length(MRA[[n_layers]]$n_dims) > 1) {
    #     Q[[n_layers]] <- do.call(bdiag.spam, Q[[n_layers]])
    # }
    # class(Q[[n_layers]]) <- "spam"

    # initialize the middle layers
    if (n_layers > 1) {
        for (j in n_layers:2) {
            # alpha_x[[j-1]] <- rnorm(ncol(MRA[[j]]$W), 0, 1)
            # alpha_y[[j-1]] <- rnorm(ncol(MRA[[j]]$W), 0, 1)
            alpha_x[[j-1]] <- rnorm(ncol(MRA[[j]]$W), 0, 0.1)
            alpha_y[[j-1]] <- rnorm(ncol(MRA[[j]]$W), 0, 0.1)

            MRA[[j-1]] <- eval_basis(cbind(MRA[[j]]$W %*% alpha_x[[j-1]], MRA[[j]]$W %*% alpha_y[[j-1]]), grid, use_spam = use_spam, ncores = ncores, nchunks = nchunks)
            Q[[j-1]] <- make_Q(sqrt(MRA[[j-1]]$n_dims),
                                      phi=rep(0.9, length(MRA[[j-1]]$n_dims)), use_spam = use_spam)

            # Q[[j-1]] <- make_Q_alpha_2d(sqrt(MRA[[j-1]]$n_dims),
            #                           phi=rep(0.9, length(MRA[[j-1]]$n_dims)))
            # if (length(MRA[[j-1]]$n_dims) > 1) {
            #     Q[[j-1]] <- do.call(bdiag.spam, Q[[j-1]])
            # }
            # class(Q[[j-1]]) <- "spam"
        }
    }

    # initialize the top layer

    # if (is.null(alpha)) {
    alpha <- rnorm(ncol(MRA[[1]]$W), 0, 0.1)
    # alpha <- rnorm(ncol(MRA[[1]]$W), 0, 1)
    # }

    z <- MRA[[1]]$W %*% alpha

    epsilon <- rnorm(N, 0, sigma)
    y <- z + epsilon
    return(list(y = y, z = z, locs = locs, grid = grid, MRA = MRA, Q = Q,
                alpha = alpha, alpha_x = alpha_x, alpha_y = alpha_y))
}
