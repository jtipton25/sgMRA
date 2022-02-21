# calc_ddiff <- function(locs, grid, MRA) {
#     # only calculate the differenxe for the non-zero rows
#     diffs <- 2*(locs[MRA$D$ind[, 1], ] - grid$locs_grid[[1]][MRA$D$ind[, 2], ])
#     MRA$D$diffsx <- diffs$x
#     ddistx <- spam(MRA$D[c("ind", "diffsx")], nrow = MRA$D$da[1], ncol = MRA$D$da[2])
#     MRA$D$diffsy <- diffs$y
#     ddisty <- spam(MRA$D[c("ind", "diffsy")], nrow = MRA$D$da[1], ncol = MRA$D$da[2])
#     return(list(ddistx=ddistx, ddisty=ddisty))
# }

update_deep_mra <- function(y, locs, grid, n_layers,
                            MRA,
                            alpha,
                            alpha_x, alpha_y,
                            Q,
                            i,
                            m, v,
                            learn_rate) {

    N <- length(y)
    # first layer
    delta <- drop(- 1/N * (y - MRA[[1]]$W %*% alpha))

    # second layer
    if (n_layers > 1) {
        delta_x[[1]] <- delta * ((MRA[[1]]$dW * MRA[[1]]$ddistx) %*% alpha_x[[1]])# check these
        delta_y[[1]] <- delta * ((MRA[[1]]$dW * MRA[[1]]$ddisty) %*% alpha_y[[1]])# check these
    }
    if (n_layers >= 2) {
        for (j in 2:n_layers) {
            delta_x[[j]] <- (delta_x[[j-1]] + delta_y[[j-1]]) * (MRA[[j]]$dW * MRA[[j]]$ddistx) %*% alpha_x[[j]]
            delta_y[[j]] <- (delta_x[[j-1]] + delta_y[[j-1]]) * (MRA[[j]]$dW * MRA[[j]]$ddistx) %*% alpha_y[[j]]
        }
    }

    # update the gradient with adam
    grad <- list(t(MRA$W) %*% delta + Q %*% alpha / N,
                 t(MRA1$W) %*% delta_x1 + Q1 %*% alpha_x1 / N,
                 t(MRA1$W) %*% delta_y1 + Q1 %*% alpha_y1 / N,
                 t(MRA2$W) %*% delta_x2 + Q2 %*% alpha_x2 / N,
                 t(MRA2$W) %*% delta_y2 + Q2 %*% alpha_y2 / N)

    adam_out <- adam(i, grad, m, v)
    alpha    <- alpha    - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
    alpha_x1 <- alpha_x1 - learn_rate * adam_out$m_hat[[2]] / (sqrt(adam_out$v_hat[[2]]) + adam_out$epsilon)
    alpha_y1 <- alpha_y1 - learn_rate * adam_out$m_hat[[3]] / (sqrt(adam_out$v_hat[[3]]) + adam_out$epsilon)
    alpha_x2 <- alpha_x2 - learn_rate * adam_out$m_hat[[4]] / (sqrt(adam_out$v_hat[[4]]) + adam_out$epsilon)
    alpha_y2 <- alpha_y2 - learn_rate * adam_out$m_hat[[5]] / (sqrt(adam_out$v_hat[[5]]) + adam_out$epsilon)


    # # decaying learning rate
    # alpha <- alpha - learn_rate * (t(MRA$W) %*% delta + Q %*% alpha / N)
    # alpha_x1 <- alpha_x1 - learn_rate * (t(MRA1$W) %*% delta_x1 + Q1 %*% alpha_x1 / N)
    # alpha_y1 <- alpha_y1 - learn_rate * (t(MRA1$W) %*% delta_y1 + Q1 %*% alpha_y1)

    # alpha <- alpha - learn_rate * (t(MRA$W) %*% delta)
    # alpha_x1 <- alpha_x1 - learn_rate * (t(MRA1$W) %*% delta_x1)
    # alpha_y1 <- alpha_y1 - learn_rate * (t(MRA1$W) %*% delta_y1)
    # # alpha_x2 <-
    #     alpha_x2 - learn_rate * t(MRA2$W) %*% delta_x2
    # # alpha_y2 <-
    #     alpha_y2 - learn_rate * t(MRA2$W) %*% delta_y2
    # # alpha_x1 <-
    #     alpha_x1 - learn_rate * t(MRA1$W) %*% delta_x1
    # # alpha_y1 <-
    #     alpha_y1 - learn_rate * t(MRA1$W) %*% delta_y1

    # the forward pass on the sgMRA
    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid)
    MRA <- eval_basis(cbind(MRA1$W %*% alpha_x1, MRA1$W %*% alpha_y1), grid)

    return(list(alpha = alpha, alpha_x1 = alpha_x1, alpha_y1 = alpha_y1,
                alpha_x2 = alpha_x2, alpha_y2 = alpha_y2,
                MRA = MRA, MRA1 = MRA1,
                m=adam_out$m, v=adam_out$v))

}


fit_sgd <- function(y, locs, grid,
                    n_layers = 3,
                    # alpha = NULL,
                    # alpha_x1 = NULL,
                    # alpha_y1 = NULL,
                    # alpha_x2 = NULL,
                    # alpha_y2 = NULL,
                    learn_rate, n_iter=500,
                    # rate_schedule = NULL,
                    n_message = 50) {

    # if (is.null(rate_schedule)) {
    #     rate_schedule <- rep(learn_rate, n_iter)
    #     # seq(0.5, 0.05, length = 5), each = ceiling(n_iter / 5)),
    # }

    N <- length(y)


    # initialize the MRA parameters here
    MRA <- vector(mode = 'list', length = n_layers)
    Q <- vector(mode = 'list', length = n_layers)
    alpha_x <- vector(mode = 'list', length = n_layers-1)
    alpha_y <- vector(mode = 'list', length = n_layers-1)

    # initialize the bottom layer
    MRA[[n_layers]] <- eval_basis(locs, grid)
    Q[[n_layers]] <- make_Q_alpha_2d(sqrt(MRA[[n_layers]]$n_dims),
                                     phi=rep(0.9, length(MRA[[n_layers]]$n_dims)))
    if (length(MRA[[n_layers]]$n_dims) > 1) {
        Q[[n_layers]] <- do.call(bdiag.spam, Q[[n_layers]])
    }
    class(Q[[n_layers]]) <- "spam"

    # initialize the middle layers
    if (n_layers > 2) {
        # if (is.null(alpha_x[[n_layers-1]])) {
        alpha_x[[n_layers-1]] <- rnorm(ncol(MR[[n_layers]]$W), 0, 0.1)
        # }
        # if (is.null(alpha_y[[n_layers-1]])) {
        alpha_y[[n_layers-1]] <- rnorm(ncol(MRA[[n_layers]]$W), 0, 0.1)
        # }

        for (j in (n_layers-1):2) {
            MRA[[j]] <- eval_basis(cbind(MRA[[j+1]]$W %*% alpha_x[[j]], MRA[[j+1]]$W %*% alpha_y[[j]]), grid)
            Q[[j]] <- make_Q_alpha_2d(sqrt(MRA[[j]]$n_dims),
                                      phi=rep(0.9, length(MRA[[j]]$n_dims)))
            if (length(MRA[[j]]$n_dims) > 1) {
                Q[[j]] <- do.call(bdiag.spam, Q[[j]])
            }
            class(Q[[j]]) <- "spam"

            # if (is.null(alpha_x[[j-1]])) {
            alpha_x[[j-1]] <- rnorm(ncol(MR[[j]]$W), 0, 0.1)
            # }
            # if (is.null(alpha_y[[j-1]])) {
            alpha_y[[j-1]] <- rnorm(ncol(MRA[[j]]$W), 0, 0.1)
            # }
        }
    }

    # initialize the top layer
    MRA[[1]] <- eval_basis(cbind(MRA[[2]]$W %*% alpha_x[[1]], MRA[[2]]$W %*% alpha_y[[1]]), grid)
    Q[[1]] <- make_Q_alpha_2d(sqrt(MRA[[1]]$n_dims),
                              phi=rep(0.9, length(MRA[[1]]$n_dims)))
    if (length(MRA[[1]]$n_dims) > 1) {
        Q[[1]] <- do.call(bdiag.spam, Q[[1]])
    }
    class(Q[[1]]) <- "spam"
    # top layer (if there is at least one deep layer)

    if (is.null(alpha)) {
        alpha <- rnorm(ncol(MRA[[1]]$W), 0, 0.1)
    }

    # initialize the loss
    loss <- rep(NA, n_iter)

    # initialize adam optimization
    m <- vector(mode='list', length = 1 + 2 * (n_layers-1)) # alpha, alpha_x1, and alpha_y1
    v <- vector(mode='list', length = 1 + 2 * (n_layers-1)) # alpha, alpha_x1, and alpha_y1
    m[[1]] <- rep(0, length(alpha))
    for (j in 1:(n_layers-1)) {
        m[[2*j + 1]] <- rep(0, length(alpha_x[[j]]))
        m[[2*j + 2]] <- rep(0, length(alpha_y[[j]]))
        v[[2*j + 1]] <- rep(0, length(alpha_x[[j]]))
        v[[2*j + 2]] <- rep(0, length(alpha_y[[j]]))
    }

    for (i in 1:n_iter) {


        pars <- update_deep_mra(y, locs, grid, n_layers,
                                MRA,
                                alpha,
                                alpha_x,
                                alpha_y,
                                Q,
                                i,
                                m, v,
                                learn_rate)
        # rate_schedule[i])
        alpha <- pars$alpha
        alpha_x <- pars$alpha_x
        alpha_y <- pars$alpha_y
        MRA <- pars$MRA
        m <- pars$m
        v <- pars$v

        loss[i] <- 1 / N * sum((y - MRA[[1]]$W %*% alpha)^2)
        if (i %% n_message == 0) {
            message("iteration i = ", i, " loss = ", loss[i])
        }
    }

    return(list(alpha = alpha, alpha_x = alpha_x, alpha_y = alpha_y,
                MRA = MRA, loss = loss))
}



