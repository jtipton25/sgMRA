make_D_ind <- function(MRA) {
    n_layers <- length(MRA$n_dims)
    D_ind <- vector(mode='list', length = n_layers)
    if (n_layers == 1) {
        D_ind <- MRA$D[[1]]$ind
    } else {
        D_ind[[1]] <- MRA$D[[1]]$ind
        for (j in 2:n_layers) {
            D_ind[[j]] <- MRA$D[[j]]$ind + cbind(rep(0, nrow(MRA$D[[j]]$ind)), sum(MRA$n_dims[1:(j-1)]))
        }
        D_ind = do.call(rbind, D_ind)
    }
    return(D_ind)
}

sparse_outer_W <- function(MRA, alpha, y) {

    N <- length(y)
    N_grid <- length(alpha)
    D_ind = make_D_ind(MRA)
    devs <- 1 / N * (MRA$W %*% alpha - y)
    delta_W_sparse = devs[D_ind[, 1]] * alpha[D_ind[, 2]]
    delta_W_sparse <- sparseMatrix(i=D_ind[, 1], j=D_ind[, 2], x=delta_W_sparse, dims=c(N, N_grid))
    return(delta_W_sparse)
}

sparse_outer_layers <- function(delta, alpha, MRA) {

    N <- length(delta)
    N_grid <- length(alpha)
    D_ind = make_D_ind(MRA)
    delta_sparse = delta[D_ind[, 1]] * alpha[D_ind[, 2]]
    delta_sparse <- sparseMatrix(i=D_ind[, 1], j=D_ind[, 2], x=delta_sparse, dims=c(N, N_grid))
    return(delta_sparse)
}

update_deep_mra <- function(y, locs, grid, MRA, MRA1, MRA2,
                            alpha,
                            alpha_x1, alpha_y1,
                            alpha_x2, alpha_y2,
                            Q,
                            Q1,
                            Q2,
                            i,
                            m, v,
                            learn_rate,
                            penalized,
                            use_spam,
                            sparse_outer,
                            noisy,
                            ncores,
                            nchunks) {

    N <- length(y)
    # first layer w.r.t alpha
    delta <- drop(1 / N * (MRA$W %*% alpha - y))
    # first layer w.r.t W
    if (sparse_outer) {
        delta_W <- sparse_outer_W(MRA, alpha, y)
    } else {
        delta_W <- (1 / N * (MRA$W %*% alpha - y)) %*% t(alpha)
    }

    # second layer
    delta_x1 <- rowSums(delta_W * (MRA$dW * MRA$ddistx))
    delta_y1 <- rowSums(delta_W * (MRA$dW * MRA$ddisty))

    # third layer
    if (sparse_outer) {
        delta2 <- sparse_outer_layers(delta_x1, alpha_x1, MRA1) +
            sparse_outer_layers(delta_y1, alpha_y1, MRA1)
    } else {
        delta2 <- delta_x1 %*% t(alpha_x1) + delta_y1 %*% t(alpha_y1) # chain rule gradient tape
    }

    delta_x2 <- rowSums(delta2 * MRA1$dW * MRA1$ddistx)
    delta_y2 <- rowSums(delta2 * MRA1$dW * MRA1$ddisty)

    # update the gradient with adam with a penalty
    if (penalized) {
        if (noisy) {
            grad <- list(t(MRA$W) %*% delta + Q %*% alpha / N + rnorm(length(alpha), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA1$W) %*% delta_x1 + Q1 %*% alpha_x1 / N + rnorm(length(alpha_x1), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA1$W) %*% delta_y1 + Q1 %*% alpha_y1 / N + rnorm(length(alpha_y1), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA2$W) %*% delta_x2 + Q2 %*% alpha_x2 / N + rnorm(length(alpha_x2), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA2$W) %*% delta_y2 + Q2 %*% alpha_y2 / N + rnorm(length(alpha_y2), 0, sqrt(0.3 / (1 + i)^(0.55))))
        } else {
            grad <- list(t(MRA$W) %*% delta + Q %*% alpha / N,
                         t(MRA1$W) %*% delta_x1 + Q1 %*% alpha_x1 / N,
                         t(MRA1$W) %*% delta_y1 + Q1 %*% alpha_y1 / N,
                         t(MRA2$W) %*% delta_x2 + Q2 %*% alpha_x2 / N,
                         t(MRA2$W) %*% delta_y2 + Q2 %*% alpha_y2 / N)
        }
    } else {
        # update the gradient without a penalty term
        if (noisy) {
            grad <- list(t(MRA$W) %*% delta + rnorm(length(alpha), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA1$W) %*% delta_x1 + rnorm(length(alpha_x1), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA1$W) %*% delta_y1 + rnorm(length(alpha_y1), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA2$W) %*% delta_x2 + rnorm(length(alpha_x2), 0, sqrt(0.3 / (1 + i)^(0.55))),
                         t(MRA2$W) %*% delta_y2 + rnorm(length(alpha_y2), 0, sqrt(0.3 / (1 + i)^(0.55))))

        } else {

            grad <- list(t(MRA$W) %*% delta,
                         t(MRA1$W) %*% delta_x1,
                         t(MRA1$W) %*% delta_y1,
                         t(MRA2$W) %*% delta_x2,
                         t(MRA2$W) %*% delta_y2)
        }
    }

    adam_out <- adam(i, grad, m, v)
    alpha    <- alpha    - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
    alpha_x1 <- alpha_x1 - learn_rate * adam_out$m_hat[[2]] / (sqrt(adam_out$v_hat[[2]]) + adam_out$epsilon)
    alpha_y1 <- alpha_y1 - learn_rate * adam_out$m_hat[[3]] / (sqrt(adam_out$v_hat[[3]]) + adam_out$epsilon)
    alpha_x2 <- alpha_x2 - learn_rate * adam_out$m_hat[[4]] / (sqrt(adam_out$v_hat[[4]]) + adam_out$epsilon)
    alpha_y2 <- alpha_y2 - learn_rate * adam_out$m_hat[[5]] / (sqrt(adam_out$v_hat[[5]]) + adam_out$epsilon)


    # the forward pass on the sgMRA
    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid, use_spam=use_spam, ncores = ncores, nchunks = nchunks)
    MRA <- eval_basis(cbind(MRA1$W %*% alpha_x1, MRA1$W %*% alpha_y1), grid, use_spam=use_spam, ncores = ncores, nchunks = nchunks)

    return(list(alpha = alpha, alpha_x1 = alpha_x1, alpha_y1 = alpha_y1,
                alpha_x2 = alpha_x2, alpha_y2 = alpha_y2,
                MRA = MRA, MRA1 = MRA1,
                m=adam_out$m, v=adam_out$v))

}


#' Title
#'
#' @param y The data
#' @param locs An N x 2 matrix of spatial locations
#' @param grid A grid object that is the output of \code{make_grid}
#' @param alpha If specified, the top layer MRA parameters
#' @param alpha_x1 If specified, the first hidden layer MRA parameters for the x-axis variable
#' @param alpha_y1 If specified, the first hidden layer MRA parameters for the y-axis variable
#' @param alpha_x2 If specified, the second hidden layer MRA parameters for the x-axis variable
#' @param alpha_y2 If specified, the second hidden layer MRA parameters for the y-axis variable
#' @param learn_rate The gradient descent learning rate
#' @param rate_schedule If specified, the gradient descent learning rate schedule in decreasing values.
#' @param n_iter The number of gradient descent iterations
#' @param n_message The number of iterations between which to output a message
#' @param penalized Fit using a penalty term
#' @param plot_during_fit Plot the current parameter states every \code{n_message} iterations
#' @param use_spam Whether to use the spam (\code{use_spam = TRUE}) or Matrix (\code{use_spam = FALSE}) package for sparse matrices
#' @param adam_pars The adam parameter state to allow restarting the model
#' @param sparse_outer If \code{TRUE}, calculate the outer product in a sparse format. For all but the smallest models, this should be TRUE. I should make this automatic going forward
#' @param noisy If \code{TRUE}, add random noise to the gradient.
#' @param ncores The number of cores to use for parallelization
#' @param nchunks The number of chunks to divide the distance calculation into. The default argument of NULL will use the same number of chunks as the number of cores.

#'
#' @return
#' @export
#'
#' @import Matrix spam
#'
fit_sgd <- function(y,
                    locs,
                    grid,
                    alpha = NULL,
                    alpha_x1 = NULL,
                    alpha_y1 = NULL,
                    alpha_x2 = NULL,
                    alpha_y2 = NULL,
                    learn_rate = 0.001,
                    rate_schedule=NULL,
                    n_iter=500,
                    n_message = 50,
                    penalized = FALSE,
                    plot_during_fit = FALSE,
                    use_spam = FALSE,
                    adam_pars = NULL,
                    sparse_outer = TRUE,
                    noisy=TRUE,
                    ncores        = 1L,
                    nchunks       = NULL) {


    if (is.null(rate_schedule)) {
        rate_schedule <- rep(learn_rate, n_iter)
    }

    N <- length(y)


    # initialize the MRA parameters here
    MRA2 <- eval_basis(locs, grid, use_spam=use_spam, ncores = ncores, nchunks = nchunks)
    Q2 <- make_Q(sqrt(MRA2$n_dims), phi=rep(0.9, length(MRA2$n_dims)), use_spam=use_spam)
    if (!use_spam) {
        CH2 <- Cholesky(Q2)
    }
    # Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=rep(0.9, length(MRA2$n_dims)), use_spam=use_spam)
    # if (length(MRA2$n_dims) > 1) {
    #     for (m in 1:M){
    #         Q2[[m]] <- 2^(2*(m-1)) * Q2[[m]]
    #     }
    #     if (use_spam) {
    #         Q2 <- do.call(bdiag.spam, Q2)
    #     } else {
    #         Q2 <- do.call(bdiag, Q2)
    #     }
    # }
    # if (use_spam) {
    #     class(Q2) <- "spam"
    # } else {
    #     class(Q2) <- "dgCMatrix"
    #     CH2 <- Cholesky(Q2)
    # }

    if (is.null(alpha_x2)) {
        # alpha_x2 <- rnorm(ncol(MRA2$W), 0, 0.1)
        # alpha_x2 <- rnorm(ncol(MRA2$W), 0, 1)
        if (use_spam) {
            alpha_x2 <- drop(rmvnorm.prec(1, rep(0, ncol(MRA2$W)), Q2)) * 0.1
        } else {
            alpha_x2 <- drop(sparseMVN::rmvn.sparse(1, rep(0, ncol(MRA2$W)), CH2, prec = TRUE)) * 0.1
        }
    }
    if (is.null(alpha_y2)) {
        # alpha_y2 <- rnorm(ncol(MRA2$W), 0, 0.1)
        # alpha_y2 <- rnorm(ncol(MRA2$W), 0, 1)
        if (use_spam) {
            alpha_y2 <- drop(rmvnorm.prec(1, rep(0, ncol(MRA2$W)), Q2)) * 0.1
        } else {
            alpha_y2 <- drop(sparseMVN::rmvn.sparse(1, rep(0, ncol(MRA2$W)), CH2, prec=TRUE)) * 0.1
        }
    }

    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid, use_spam=use_spam, ncores = ncores, nchunks = nchunks)
    Q1 <- make_Q(sqrt(MRA1$n_dims), phi=rep(0.9, length(MRA1$n_dims)), use_spam=use_spam)
    if (!use_spam) {
        CH1 <- Cholesky(Q1)
    }
    # Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=rep(0.9, length(MRA1$n_dims)), use_spam=use_spam)
    # if (length(MRA1$n_dims) > 1) {
    #     for (m in 1:M){
    #         Q1[[m]] <- 2^(2*(m-1)) * Q1[[m]]
    #     }
    #     if (use_spam) {
    #         Q1 <- do.call(bdiag.spam, Q1)
    #     } else {
    #         Q1 <- do.call(bdiag, Q1)
    #     }
    #
    # }
    # if (use_spam) {
    #     class(Q1) <- "spam"
    # } else {
    #     class(Q1) <- "dgCMatrix"
    #     CH1 <- Cholesky(Q1)
    # }

    if (is.null(alpha_x1)) {
        # alpha_x1 <- rnorm(ncol(MRA1$W), 0, 0.1)
        # alpha_x1 <- rnorm(ncol(MRA1$W), 0, 1)
        if (use_spam) {
            alpha_x1 <- drop(rmvnorm.prec(1, rep(0, ncol(MRA1$W)), Q1)) * 0.1
        } else {
            alpha_x1 <- drop(sparseMVN::rmvn.sparse(1, rep(0, ncol(MRA1$W)), CH1, prec = TRUE)) * 0.1
        }
    }
    if (is.null(alpha_y1)) {
        # alpha_y1 <- rnorm(ncol(MRA1$W), 0, 0.1)
        # alpha_y1 <- rnorm(ncol(MRA1$W), 0, 1)
        if (use_spam) {
            alpha_y1 <- drop(rmvnorm.prec(1, rep(0, ncol(MRA1$W)), Q1)) * 0.1
        } else {
            alpha_y1 <- drop(sparseMVN::rmvn.sparse(1, rep(0, ncol(MRA1$W)), CH1, prec=TRUE)) * 0.1
        }
    }

    MRA <- eval_basis(cbind(MRA1$W %*% alpha_x1, MRA1$W %*% alpha_y1), grid, use_spam=use_spam, ncores = ncores, nchunks = nchunks)
    Q <- make_Q(sqrt(MRA$n_dims), phi=rep(0.9, length(MRA$n_dims)), use_spam=use_spam)
    if (!use_spam) {
        CH <- Cholesky(Q)
    }
    # Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=rep(0.9, length(MRA1$n_dims)), use_spam=use_spam)
    # if (length(MRA$n_dims) > 1) {
    #     for (m in 1:M){
    #         Q[[m]] <- 2^(2*(m-1)) * Q[[m]]
    #     }
    #     if (use_spam) {
    #         Q <- do.call(bdiag.spam, Q)
    #     } else {
    #         Q <- do.call(bdiag, Q)
    #     }
    # }
    #
    # if (use_spam) {
    #     class(Q) <- "spam"
    # } else {
    #     class(Q) <- "dgCMatrix"
    #     CH <- Cholesky(Q)
    # }

    if (is.null(alpha)) {
        # alpha <- rnorm(ncol(MRA$W), 0, 0.1)
        # alpha <- rnorm(ncol(MRA$W), 0, 1)
        if (use_spam) {
            alpha <- drop(rmvnorm.prec(1, rep(0, ncol(MRA$W)), Q)) * 0.1
        } else {
            alpha <- drop(sparseMVN::rmvn.sparse(1, rep(0, ncol(MRA$W)), CH, prec=TRUE)) * 0.1
        }
    }

    # initialize the loss
    loss <- rep(NA, n_iter)

    message("Initializing the model, the initialization loss is = ", 1 / (2 * N) * sum((y - MRA$W %*% alpha)^2))

    # plot the data
    if (plot_during_fit) {
        dat <- data.frame(x = locs$x, y = locs$y, z = drop(z), y_obs = drop(y_obs))
        p1 <- ggplot(dat, aes(x = x, y = y, fill = z)) +
            geom_raster() +
            scale_fill_viridis_c()
        dat <- data.frame(x = locs$x, y = locs$y,
                          layer = rep(c(1, 1, 2, 2, 3), each=N),
                          group = rep(c("x", "y", "x", "y", "z"), each = N),
                          z = c(drop(MRA2$W %*% alpha_x2), drop(MRA2$W %*% alpha_y2),
                                drop(MRA1$W %*% alpha_x1), drop(MRA1$W %*% alpha_y1),
                                drop(MRA$W %*% alpha)))
        p_layers_fit <- ggplot(dat, aes(x, y, fill=z)) +
            geom_raster() +
            scale_fill_viridis_c() +
            facet_grid(layer ~ group) +
            ggtitle("initialized fitted layers")
        print(p_layers_fit / p1)
    }
    message("Fitting the model for ", n_iter, " iterations")

    # initialize adam optimization
    # m <- vector(mode='list', length = 3) # alpha, alpha_x1, and alpha_y1
    # v <- vector(mode='list', length = 3) # alpha, alpha_x1, and alpha_y1
    if (!is.null(adam_pars)) {
        m <- adam_pars$m
        v <- adam_pars$v
    } else {
        m <- vector(mode='list', length = 5)
        v <- vector(mode='list', length = 5)
        m[[1]] <- rep(0, length(alpha))
        m[[2]] <- rep(0, length(alpha_x1))
        m[[3]] <- rep(0, length(alpha_y1))
        m[[4]] <- rep(0, length(alpha_x2))
        m[[5]] <- rep(0, length(alpha_y2))
        v[[1]] <- rep(0.01, length(alpha))
        v[[2]] <- rep(0.01, length(alpha_x1))
        v[[3]] <- rep(0.01, length(alpha_y1))
        v[[4]] <- rep(0.01, length(alpha_x2))
        v[[5]] <- rep(0.01, length(alpha_y2))
    }

    for (i in 1:n_iter) {


        pars <- update_deep_mra(y = y,
                                locs = locs,
                                grid = grid,
                                MRA = MRA,
                                MRA1 = MRA1,
                                MRA2 = MRA2,
                                alpha,
                                alpha_x1, alpha_y1,
                                alpha_x2, alpha_y2,
                                Q, Q1, Q2,
                                i,
                                m, v,
                                learn_rate   =  rate_schedule[i],
                                penalized    = penalized,
                                use_spam     = use_spam,
                                sparse_outer = sparse_outer,
                                noisy        = noisy,
                                ncores       = ncores,
                                nchunks      = nchunks)
        alpha <- pars$alpha
        alpha_x1 <- pars$alpha_x1
        alpha_y1 <- pars$alpha_y1
        alpha_x2 <- pars$alpha_x2
        alpha_y2 <- pars$alpha_y2
        MRA <- pars$MRA
        MRA1 <- pars$MRA1
        m <- pars$m
        v <- pars$v

        loss[i] <- 1 / (2 * N) * sum((y - MRA$W %*% alpha)^2)
        if (i %% n_message == 0) {
            message("iteration i = ", i, " loss = ", loss[i])
            if (plot_during_fit) {
                # examine the fitted layers
                dat <- data.frame(x = locs$x, y = locs$y,
                                  layer = rep(c(1, 1, 2, 2, 3), each=N),
                                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                                  z = c(drop(MRA2$W %*% alpha_x2), drop(MRA2$W %*% alpha_y2),
                                        drop(MRA1$W %*% alpha_x1), drop(MRA1$W %*% alpha_y1),
                                        drop(MRA$W %*% alpha)))
                p_layers_fit <- ggplot(dat, aes(x, y, fill=z)) +
                    geom_raster() +
                    scale_fill_viridis_c() +
                    facet_grid(layer ~ group) +
                    ggtitle(paste0("fitted layers, iteration ", i))
                print(p_layers_fit / p1)
            }

        }
    }

    return(list(alpha = alpha, alpha_x1 = alpha_x1, alpha_y1 = alpha_y1,
                alpha_x2 = alpha_x2, alpha_y2 = alpha_y2,
                MRA = MRA, MRA1 = MRA1, MRA2 = MRA2, loss = loss,
                adam_pars = list(m = pars$m, v = pars$v)))
}



