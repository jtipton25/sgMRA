#' Make the deep MRA grid
#'
#' @param locs An N x 2 matrix of spatial locations
#' @param min_x The minimum value of the MRA grid in the x axis
#' @param min_y The minimum value of the MRA grid in the y axis
#' @param max_x The maximum value of the MRA grid in the x axis
#' @param max_y The maximum value of the MRA grid in the y axis
#' @param M The number of resolutions.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param n_padding The number of additional boundary points to add on each boundary. For example, n_padding = 5 will add 5 boundary knots to the both the left  and right side of the grid).
#' @param n_neighbors The expected number of neighbors for each interior basis function. This determines the basis radius parameter.
#' @param max_points The maximum number of points in the
#' @param basis_type The basis type. Currently only the "wendland" basis is supported
#' @param use_spam Whether to use spam or Matrix for sparse matrix
#'
#' @return A grid object with \code{locs_grid} for the MRA grid and \code{radius} for the thresholded distance radius
#' @export
#'
#'
make_grid <- function(locs,
    min_x = -6,
    min_y = -6,
    max_x = 6,
    max_y = 6,
    M             = 4,
    n_coarse_grid = 10,
    n_padding     = 5L,
    n_neighbors   = 68,
    max_points    = NULL,
    basis_type    = "wendland"
) {


    N <- nrow(locs)

    ## finest grid is the smaller of n_max_fine_grid
    ## or the largest power of 2 larger than N
    n_grid <- n_coarse_grid * 2^(1:M - 1)
    if (min(n_grid) < 4) {
        stop("There are too many resolutions to form a reliable grid. Reduce M and try again.")
    }

    ## define radius so that each basis function covers approximately n_neighbors
    area_domain      <- (max_x - min_x) * (max_y - min_y)
    density_domain   <- max(n_grid^2) / area_domain
    radius_fine_grid <- sqrt(n_neighbors / (density_domain * base::pi))
    radius           <- radius_fine_grid * (2^(M:1 - 1))

    ## generate a set of grid locations for the basis
    locs_grid  <- vector(mode = "list", length = M)
    # W            <- vector(mode = "list", length = M)
    delta_x    <- vector(mode = "list", length = M)
    delta_y    <- vector(mode = "list", length = M)


    for (m in 1:M) {

        ## right now assuming a 2D grid -- can generalize later
        seq_x <- seq(
            min_x,
            max_x,
            length.out = n_grid[m]
        )
        seq_y <- seq(
            min_y,
            max_y,
            length.out = n_grid[m]
        )
        delta_x[[m]] <- seq_x[2] - seq_x[1]
        delta_y[[m]] <- seq_y[2] - seq_y[1]
        seq_x <- c(
            min(seq_x) - delta_x[[m]] * (n_padding:1),
            seq_x,
            max(seq_x) + delta_x[[m]] * (1:n_padding)
        )
        seq_y <- c(
            min(seq_y) - delta_y[[m]] * (n_padding:1),
            seq_y,
            max(seq_y) + delta_y[[m]] * (1:n_padding)
        )

        locs_grid[[m]] <- expand.grid(seq_x, seq_y)

    }

    return(list(locs_grid = locs_grid, radius = radius,
                delta_x = delta_x, delta_y = delta_y))

}



#' Evaluate the MRA basis
#'
#' @param locs An N x 2 matrix of spatial locations
#' @param grid A grid object that is the output of \code{make_grid}
#' @param basis_type The basis function type. Currently only "wendland" is accepted
#' @param use_spam Whether to use the spam (\code{use_spam = TRUE}) or Matrix (\code{use_spam = FALSE}) package for sparse matrices
#' @param ncores The number of cores to use for parallelization
#' @param nchunks The number of chunks to divide the distance calculation into. The default argument of NULL will use the same number of chunks as the number of cores.
#'
#' @return
#' @export
#'
#' @importFrom BayesMRA make_basis
#'
eval_basis <- function(
    locs,
    grid,
    basis_type    = "wendland",
    use_spam      = TRUE,
    ncores        = 1L,
    nchunks       = NULL
) {
    ##
    ## check inputs
    ##

    # if (is.null(nrow(locs))) {
    #     stop("locs must be a numeric matrix with N rows and 2 columns")
    # }
    # N <- nrow(locs)
    #
    # if (!is_numeric_matrix(locs, N, 2)) {
    #     stop("locs must be a numeric matrix with N rows and 2 columns")
    # }
    # if (!is_positive_integer(M, 1)) {
    #     stop("the number of resolutions M must be a positive integer")
    # }
    # if (!is_positive_integer(n_coarse_grid, 1)) {
    #     stop("n_coarse_grid must be a positive integer")
    # }
    # if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
    #     stop("use_spam must be either TRUE or FALSE")
    # }
    # if (use_spam == FALSE) {
    #     stop("The Matrix package is not currently supported")
    # }
    # if (!(basis_type %in% c("wendland"))) {
    #     stop('The only basis function type allowed is "wendland"')
    # }

    N <- nrow(locs)
    ncores <- as.integer(ncores)
    if (!is.null(nchunks)) {
        nchunks <- as.integer(nchunks)
    }

    ## Assign as many gridpoints (approximately) as data
    # n_grid <- ceiling(sqrt(N / 2^(M:1 - 1)))

    ## finest grid is the smaller of n_max_fine_grid
    ## or the largest power of 2 larger than N
    # n_grid <- ceiling(min(n_max_fine_grid, 2^ceiling(log(N, base = 2)))^(0.5) / 2^(M:1 - 1))

    M <- length(grid$locs_grid)
    D            <- vector(mode = "list", length = M)
    W            <- vector(mode = "list", length = M)
    dW           <- vector(mode = "list", length = M)
    ddistx       <- vector(mode = "list", length = M)
    ddisty       <- vector(mode = "list", length = M)


    for (m in 1:M) {

        # rewriten to include dW and ddist functions for more efficiency,
        # TODO: rewrite to include basis and dbasis -- not worth it as this doesn't add much speed


        if (ncores == 1) {
            D[[m]] <- distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                                          radius = grid$radius[m])
        } else {
            # TODO: make the indexing/renaming more efficient
            D_tmp <- do.call(rbind,
                              distance_near_chunk_cpp(as.matrix(locs),
                                                      as.matrix(grid$locs_grid[[m]]),
                                                      radius = grid$radius[m],
                                                      nchunks = nchunks,
                                                      ncores = ncores))
            # TODO: make the indexing/renaming more efficient
            D[[m]]$ind <- as.matrix(D_tmp[, c(1, 2)])
            D[[m]]$V <- D_tmp[, 3]
            D[[m]]$ddistx <- D_tmp[, 4]
            D[[m]]$ddisty <- D_tmp[, 5]
        }

        D[[m]]$basis <- make_basis(D[[m]]$V, grid$radius[m], basis_type='wendland')
        D[[m]]$dbasis <- dwendland_basis(D[[m]]$V, grid$radius[m])

        N_grid <- nrow(grid$locs_grid[[m]])

        if (use_spam) {
            ## use the spam sparse matrix package
            W[[m]] <- spam(D[[m]][c('ind', 'basis')], nrow =N, ncol = N_grid)
            dW[[m]] <- spam(D[[m]][c("ind", "dbasis")], nrow = N, ncol = N_grid)
            ddistx[[m]] <- spam(D[[m]][c("ind", "ddistx")], nrow = N, ncol = N_grid)
            ddisty[[m]] <- spam(D[[m]][c("ind", "ddisty")], nrow = N, ncol = N_grid)

            # }
        } else {
            # stop("The Matrix package is not currently supported")
            W[[m]] <- sparseMatrix(i = D[[m]]$ind[, 1], j=D[[m]]$ind[, 2], x = as.numeric(D[[m]]$basis), dims = c(N, N_grid))
            dW[[m]] <- sparseMatrix(i = D[[m]]$ind[, 1], j=D[[m]]$ind[, 2], x = as.numeric(D[[m]]$dbasis), dims = c(N, N_grid))
            ddistx[[m]] <- sparseMatrix(i = D[[m]]$ind[, 1], j=D[[m]]$ind[, 2], x = as.numeric(D[[m]]$ddistx), dims = c(N, N_grid))
            ddisty[[m]] <- sparseMatrix(i = D[[m]]$ind[, 1], j=D[[m]]$ind[, 2], x = as.numeric(D[[m]]$ddisty), dims = c(N, N_grid))
            # }
        }
    }

    ## calculate the basis dimensions
    ## n_dims is the number of columns in the basis matrix for each resolution
    ## dims_idx is a vector of which
    n_dims   <- rep(NA, length(W))
    dims_idx <- c()
    for (i in 1:M) {
        n_dims[i] <- ncol(W[[i]])
        dims_idx  <- c(dims_idx, rep(i, n_dims[i]))
    }

    ## flatten the list of basis functions W to a single matrix
    # D <- do.call(cbind, D)
    W <- do.call(cbind, W)
    dW <- do.call(cbind, dW)
    ddistx <- do.call(cbind, ddistx)
    ddisty <- do.call(cbind, ddisty)

    out <- list(
        locs          = locs,
        locs_grid     = grid$locs_grid,
        W             = W,
        dW            = dW,
        ddistx        = ddistx,
        ddisty        = ddisty,
        D             = D,
        radius        = grid$radius,
        M             = M,
        n_dims        = n_dims,
        dims_idx      = dims_idx,
        n_coarse_grid = n_coarse_grid,
        use_spam      = use_spam
    )

    return(out)
}
