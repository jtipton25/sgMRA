# calculate the chain rule derivatives
make_D_nonstationary <- function(sq_devs, a, b) {
    N <- dim(sq_devs)[1]
    N_grid <- dim(sq_devs)[2]
    D <- matrix(0, N, N_grid)
    for (i in 1:N) {
        D[i, ] <- sqrt(rowSums(sq_devs[i, , ] * cbind(rep(a[i], N_grid), rep(b[i], N_grid))))
    }

    return(D)
}


d_kernel_link <- function(a, trunc) {
    out <- rep(1, length(a))
    out[a <= trunc] <- 0
    return(out)
}


# make_D_ind <- function(MRA) {
#     n_layers <- length(MRA$n_dims)
#     D_ind <- vector(mode='list', length = n_layers)
#     if (n_layers == 1) {
#         D_ind <- MRA$D[[1]]$ind
#     } else {
#         D_ind[[1]] <- MRA$D[[1]]$ind
#         for (j in 2:n_layers) {
#             D_ind[[j]] <- MRA$D[[j]]$ind + cbind(rep(0, nrow(MRA$D[[j]]$ind)), sum(MRA$n_dims[1:(j-1)]))
#         }
#         D_ind = do.call(rbind, D_ind)
#     }
#     return(D_ind)
# }

sparse_outer_W <- function(D_ind, W, alpha, y) {

    N <- length(y)
    N_grid <- length(alpha)
    devs <- 1 / N * (W %*% alpha - y)
    delta_W_sparse = devs[D_ind[, 1]] * alpha[D_ind[, 2]]
    delta_W_sparse <- sparseMatrix(i = D_ind[, 1], j = D_ind[, 2], x = delta_W_sparse, dims = c(N, N_grid))
    return(delta_W_sparse)
}


update_nonstationary_MRA <- function(y, locs, grid, sq_devs, D, radius, W, W_kernels, trunc,
                                     alpha, grid_kernels, alpha_a, alpha_b,
                                     i, m, v, learn_rate) {

    N <- length(y)
    # D_ind <- which(D < radius, arr.ind = TRUE)
    D_ind <- cbind(D@i, D@j) + 1
    delta <- drop(1 / N * (W %*% alpha - y))

    # eventually make this sparse and only keep the terms needed
    # delta_W <- (1 / N * (W %*% alpha - y)) %*% t(alpha)
    # sparse_outer_W
    delta_W <- sparse_outer_W(D_ind, W, alpha, y)

    dW <- D
    dW@x <- dwendland_basis(D@x, radius)
    # dW <- Matrix(dwendland_basis(D, radius), sparse = TRUE)
    D_inverted <- D
    D_inverted@x <- 1 / D_inverted@x

    # remove the zero distances when obseration locations are exactly on gridcell locations
    zeros_idx <- which(D@x == 0)
    D_inverted@x[zeros_idx] <- 0
    # ddista <- sq_devs[, , 1] / D
    # ddistb <- sq_devs[, , 2] / D
    ddista <- sq_devs[, , 1] * D_inverted
    ddistb <- sq_devs[, , 2] * D_inverted

    # don't need this for sparse matrix
    # # rewrite the derivative when locations are at the grid points
    # # need to think about whether this is true
    # idx_D <- which(D == 0)
    # ddista[idx_D] <- 0
    # ddistb[idx_D] <- 0


    d_tobit_a <- d_kernel_link(as.vector(W_kernels %*% alpha_a), trunc)
    d_tobit_b <- d_kernel_link(as.vector(W_kernels %*% alpha_b), trunc)

    delta_W_dW <- delta_W * dW

    delta_a <- rowSums(delta_W_dW * ddista * d_tobit_a)
    delta_b <- rowSums(delta_W_dW * ddistb * d_tobit_b)

    grad <- list(t(W) %*% delta,
                 t(W_kernels) %*% delta_a,
                 t(W_kernels) %*% delta_b)

    adam_out <- adam(i, grad, m, v)

    alpha   <- alpha   - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
    alpha_a <- alpha_a - learn_rate * adam_out$m_hat[[2]] / (sqrt(adam_out$v_hat[[2]]) + adam_out$epsilon)
    alpha_b <- alpha_b - learn_rate * adam_out$m_hat[[3]] / (sqrt(adam_out$v_hat[[3]]) + adam_out$epsilon)

    # D <- make_D_nonstationary(sq_devs, pmax(W_kernels %*% alpha_a, trunc), pmax(W_kernels %*% alpha_b, trunc))
    # eventually make this sparse
    # D <- distance_nonstationary(sq_devs,
    #                                   as.vector(pmax(W_kernels %*% alpha_a, trunc)),
    #                                   as.vector(pmax(W_kernels %*% alpha_b, trunc)))
    DD <- distance_near_nonstationary(sq_devs,
                                      as.vector(pmax(W_kernels %*% alpha_a, trunc)),
                                      as.vector(pmax(W_kernels %*% alpha_b, trunc)),
                                     radius = radius)
    # truncate to a sparse matrix
    # D[D >= radius] <- 0
    # D <- as(D, "dgTMatrix")
    D <- sparseMatrix(i = DD$ind[, 1], j = DD$ind[, 2], x = as.vector(DD$V), dims = dim(sq_devs)[1:2], repr = "T")


        # D <- distance_near_nonstationary(sq_devs,
    #                                   as.vector(pmax(W_kernels %*% alpha_a, trunc)),
    #                                   as.vector(pmax(W_kernels %*% alpha_b, trunc)),
    #                                   radius = radius)

    # W1 <- sparseMatrix(i = D@i, p = D@p, x = BayesMRA::wendland_basis(D@x, radius), dims = D@Dim)


    W <- D
    W@x <- BayesMRA::wendland_basis(D@x, radius)
    # W <- Matrix(BayesMRA::wendland_basis(D, radius), sparse = TRUE)

    return(list(alpha = alpha, alpha_a = alpha_a, alpha_b = alpha_b,
                D = D, W = W, m=adam_out$m, v=adam_out$v))

}



#' Title
#'
#' @param y The data
#' @param locs An N x 2 matrix of spatial locations
#' @param grid A grid object
#' @param grid_kernels A grid object for the spatially-varying kernels
#' @param trunc The truncation parameter for the kernel parameters
#' @param alpha If specified, the top layer MRA parameters
#' @param alpha_a If specified, the parameters for the x co-ordinate of the kernels
#' @param alpha_b If specified, the parameters for the y co-ordinate of the kernels
#' @param learn_rate The gradient descent learning rate
#' @param rate_schedule If specified, the gradient descent learning rate schedule in decreasing values.
#' @param n_iter The number of gradient descent iterations
#' @param n_message The number of iterations between which to output a message
#' @param plot_during_fit Plot the current parameter states every \code{n_message} iterations
#' @param adam_pars The adam parameter state to allow restarting the model
#'
#' @return
#' @export
#'
#' @import Matrix spam
#'
fit_nonstationary_MRA <- function(y,
                                  locs,
                                  grid,
                                  grid_kernels,
                                  radius = 0.5,
                                  trunc = 0.01,
                                  alpha = NULL,
                                  alpha_a = NULL,
                                  alpha_b = NULL,
                                  learn_rate = 0.001,
                                  rate_schedule=NULL,
                                  n_iter=500,
                                  n_message = 50,
                                  plot_during_fit = FALSE,
                                  adam_pars = NULL) {

    library(Matrix)


    if (is.null(rate_schedule)) {
        rate_schedule <- rep(learn_rate, n_iter)
    }

    N <- length(y)

    # initialize the kernel process
    D_kernels <- matrix(0, N, N_grid_kernels)
    for (i in 1:N) {
        D_kernels[i, ] <- sqrt(rowSums(sweep(grid_kernels, 2, FUN='-', locs[i, ])^2))
    }

    W_kernels <- Matrix(BayesMRA::wendland_basis(D_kernels, radius), sparse = TRUE)
    if (is.null(alpha_a)) {
        # add in penalty parameter later
        alpha_a <- rnorm(N_grid_kernels, 0, 0.1)
    }
    if (is.null(alpha_b)) {
        # add in penalty parameter later
        alpha_b <- rnorm(N_grid_kernels, 0, 0.1)
    }
    a <- pmax(drop(W_kernels %*% alpha_a), trunc)
    b <- pmax(drop(W_kernels %*% alpha_b), trunc)

    # initialize the spatial process
    sq_devs <- array(0, dim=c(N, N_grid, 2))
    for (i in 1:N) {
        # can make this parallelizable and/or Rcpp
        sq_devs[i, , ] <- sweep(grid, 2, FUN='-', locs[i, ])^2
    }

    # D <- matrix(0, N, N_grid)
    # for (i in 1:N) {
    #     D[i, ] <- sqrt(rowSums(sq_devs[i, , ] * cbind(rep(a[i], N_grid), rep(b[i], N_grid))))
    # }
    DD <- distance_near_nonstationary(sq_devs,
                                      as.vector(pmax(W_kernels %*% alpha_a, trunc)),
                                      as.vector(pmax(W_kernels %*% alpha_b, trunc)),
                                      radius = radius)
    # truncate to a sparse matrix
    D <- sparseMatrix(i = DD$ind[, 1], j = DD$ind[, 2], x = as.vector(DD$V), dims = dim(sq_devs)[1:2], repr = "T")


    # W <- Matrix(BayesMRA::wendland_basis(D, radius), sparse = TRUE)
    # sparse kernel matrix
    W <- D
    W@x <- BayesMRA::wendland_basis(D@x, radius)

    if (is.null(alpha)) {
        # add in penalty parameter later
        alpha <- rnorm(N_grid, 0, 0.1)
    }

    # initialize the loss
    loss <- rep(NA, n_iter)

    message("Initializing the model, the initialization loss is = ", 1 / (2 * N) * sum((y - W %*% alpha)^2))
    # plot the data
    if (plot_during_fit) {
        dat <- data.frame(x = locs[, 1], y = locs[, 2], z = y)
        p1 <- ggplot(dat, aes(x = x, y = y, fill = z)) +
            geom_raster() +
            scale_fill_viridis_c()
        dat <- data.frame(x = locs[, 1], y = locs[, 2], z = as.vector(W %*% alpha))
        p_fit <- ggplot(dat, aes(x, y, fill=z)) +
            geom_raster() +
            scale_fill_viridis_c() +
            ggtitle("initialized fitted process")
        print(p_fit / p1)
    }

    # loop
    message("Fitting the model for ", n_iter, " iterations")
    if (!is.null(adam_pars)) {
        m <- adam_pars$m
        v <- adam_pars$v
    } else {
        m <- vector(mode='list', length = 3)
        v <- vector(mode='list', length = 3)
        m[[1]] <- rep(0, length(alpha))
        m[[2]] <- rep(0, length(alpha_a))
        m[[3]] <- rep(0, length(alpha_b))
        v[[1]] <- rep(0.01, length(alpha))
        v[[2]] <- rep(0.01, length(alpha_a))
        v[[3]] <- rep(0.01, length(alpha_b))

    }
    # Start the gradient descent loop
    for (i in 1:n_iter) {

        pars <- update_nonstationary_MRA(y, locs, grid, sq_devs, D, radius, W, W_kernels,
                                         trunc, alpha, grid_kernels, alpha_a, alpha_b,
                                         i, m, v, rate_schedule[i])

        alpha <- pars$alpha
        alpha_a <- pars$alpha_a
        alpha_b <- pars$alpha_b
        D <- pars$D
        W <- pars$W
        m <- pars$m
        v <- pars$v

        loss[i] <- 1 / (2 * N) * sum((y - W %*% alpha)^2)
        if (i %% n_message == 0) {
            message("iteration i = ", i, " loss = ", loss[i])
            if (plot_during_fit) {
                # examine the fitted process
                dat <- data.frame(x = locs[, 1], y = locs[, 2], z = as.vector(W %*% alpha))
                p_fit <- ggplot(dat, aes(x, y, fill=z)) +
                    geom_raster() +
                    scale_fill_viridis_c() +
                    ggtitle(paste0("fitted layers, iteration ", i))
                print(p_fit / p1)
            }

        }
    }

    return(list(alpha = alpha, alpha_a = alpha_a, alpha_b = alpha_b,
                D = D, W = W, W_kernels = W_kernels, loss = loss,
                adam_pars = list(m = pars$m, v = pars$v)))
}



