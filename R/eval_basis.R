make_grid <- function(locs,
    min_x = -6, min_y = -6, max_x = 6, max_y = 6,
    M             = 4,
    n_coarse_grid = 10,
    n_padding     = 5L,
    n_neighbors   = 68,
    max_points    = NULL,
    # n_max_fine_grid = 2^12,
    # radius      = 25,
    basis_type    = "wendland",
    use_spam      = TRUE
) {


    N <- nrow(locs)

    ## Define max_points parameter
    if (is.null(max_points)) {
        max_points <- 2* N * n_neighbors
    }
    ## Assign as many gridpoints (approximately) as data
    # n_grid <- ceiling(sqrt(N / 2^(M:1 - 1)))

    ## finest grid is the smaller of n_max_fine_grid
    ## or the largest power of 2 larger than N
    # n_grid <- ceiling(min(n_max_fine_grid, 2^ceiling(log(N, base = 2)))^(0.5) / 2^(M:1 - 1))
    n_grid <- n_coarse_grid * 2^(1:M - 1)
    if (min(n_grid) < 4) {
        stop("There are too many resolutions to form a reliable grid. Reduce M and try again.")
    }

    ## define radius so that each basis function covers approximately 5 neighbors
    area_domain      <- (max_x - min_x) * (max_y - min_y)
    density_domain   <- max(n_grid^2) / area_domain
    radius_fine_grid <- sqrt(n_neighbors / (density_domain * base::pi))
    # radius           <- radius_fine_grid * (2^(1:M - 1))^2
    radius           <- radius_fine_grid * (2^(M:1 - 1))
    # radius <- radius / 2^(1:M - 1)

    ## generate a set of grid locations for the basis
    locs_grid    <- vector(mode = "list", length = M)
    out          <- vector(mode = "list", length = M)
    W            <- vector(mode = "list", length = M)

    # guess the max_points variable

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
        delta_x <- seq_x[2] - seq_x[1]
        delta_y <- seq_y[2] - seq_y[1]
        seq_x <- c(
            min(seq_x) - delta_x * (n_padding:1),
            seq_x,
            max(seq_x) + delta_x * (1:n_padding)
        )
        seq_y <- c(
            min(seq_y) - delta_y * (n_padding:1),
            seq_y,
            max(seq_y) + delta_y * (1:n_padding)
        )

        locs_grid[[m]] <- expand.grid(seq_x, seq_y)

    }

    return(list(locs_grid=locs_grid, radius = radius))

}



eval_basis <- function(
    locs,
    grid,
    n_padding     = 5L,
    n_neighbors   = 68,
    max_points    = NULL,
    # n_max_fine_grid = 2^12,
    # radius      = 25,
    basis_type    = "wendland",
    use_spam      = TRUE
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
    # if (!is_positive_integer(n_neighbors, 1)) {
    #     stop("n_neighbors must be a positive integer")
    # }
    # if (!is_positive_integer(n_coarse_grid, 1)) {
    #     stop("n_coarse_grid must be a positive integer")
    # }
    # if (!is_positive_integer(n_padding, 1)) {
    #     stop("n_padding must be a positive integer")
    # }
    # if (!(is.null(max_points) | is_positive_integer(max_points, 1))) {
    #     stop("max_points must be either NULL or a positive numeric integer")
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

    ## Define max_points parameter
    if (is.null(max_points)) {
        max_points <- 2* N * n_neighbors
    }
    ## Assign as many gridpoints (approximately) as data
    # n_grid <- ceiling(sqrt(N / 2^(M:1 - 1)))

    ## finest grid is the smaller of n_max_fine_grid
    ## or the largest power of 2 larger than N
    # n_grid <- ceiling(min(n_max_fine_grid, 2^ceiling(log(N, base = 2)))^(0.5) / 2^(M:1 - 1))

    M <- length(grid$locs_grid)
    W            <- vector(mode = "list", length = M)
    dW           <- vector(mode = "list", length = M)
    ddistx       <- vector(mode = "list", length = M)
    ddisty       <- vector(mode = "list", length = M)

    # guess the max_points variable

    for (m in 1:M) {


        # rewrite this to include dW and ddist functions for more efficiency
        # D <- fields::fields.rdist.near(locs,
        #                                grid$locs_grid[[m]],
        #                                delta = grid$radius[m],
        #                                max.points = max_points)
        #
        # D$basis <- make_basis(D$ra, grid$radius[m], basis_type = "wendland")

        D <- distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                                          radius = grid$radius[m],
                                          n_neighbors = n_neighbors)
        D$basis <- make_basis(D$V, grid$radius[m], basis_type='wendland')
        D$dbasis <- dwendland_basis(D$V, grid$radius[m])

        # if (calc_derivative) {

            # Note: this will likely fail for multiple resolutions

            # diffs <- 2*(locs[D$ind[, 1], ] - grid$locs_grid[[m]][D$ind[, 2], ]) # i think this is wrong because it doesn't include the square root
            # D$diffsx <- diffs[[1]]
            # D$diffsy <- -diffs[[2]]
            # see https://www.wolframalpha.com/input?i=d%2Fdx+sqrt%28x%5E2+%2B++y%5E2%29
            # diffs <- locs[D$ind[, 1], ] - grid$locs_grid[[m]][D$ind[, 2], ]
            # D$diffsx <- diffs[[1]] / D$ra
            # D$diffsy <- diffs[[2]] / D$ra
        # }

        if (use_spam) {
            ## use the spam sparse matrix package
            # W[[m]] <- spam(c(wendland_basis(D, grid$radius[m])), nrow = nrow(D), ncol = ncol(D))
            # W[[m]] <- spam(D[c("ind", "basis")], nrow = D$da[1], ncol = D$da[2])
            W[[m]] <- spam(D[c('ind', 'basis')], nrow =N, ncol = nrow(grid$locs_grid[[m]]))
            # if (calc_derivative) {
            dW[[m]] <- spam(D[c("ind", "dbasis")], nrow = N, ncol = nrow(grid$locs_grid[[m]]))
            ddistx[[m]] <- spam(D[c("ind", "ddistx")], nrow = N, ncol = nrow(grid$locs_grid[[m]]))
            ddisty[[m]] <- spam(D[c("ind", "ddisty")], nrow = N, ncol = nrow(grid$locs_grid[[m]]))

            # }
        } else {
            stop("The Matrix package is not currently supported")
            ## use the Matrix sparse matrix package
            W[[m]] <- Matrix(wendland_basis(D, grid$radius[m]), sparse = TRUE)
            # if (calc_derivative) {
            # note: this code is not correct for Matrix package
            dW[[m]] <- Matrix(dwendland_basis(D, grid$radius[m]), sparse = TRUE)
            ddistx[[m]] <- Matrix(D[c("ind", "ddistx")], nrow = D$da[1], ncol = nrow(grid$locs_grid[[m]]), sparse = TRUE)
            ddisty[[m]] <- Matrix(D[c("ind", "ddisty")], nrow = D$da[1], ncol = nrow(grid$locs_grid[[m]]), sparse = TRUE)

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
        n_neighbors   = n_neighbors,
        n_coarse_grid = n_coarse_grid,
        n_padding     = n_padding,
        use_spam      = use_spam
    )

    return(out)
}
