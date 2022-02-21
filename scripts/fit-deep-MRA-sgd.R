# calc_ddiff <- function(locs, grid, MRA) {
#     # only calculate the differenxe for the non-zero rows
#     diffs <- 2*(locs[MRA$D$ind[, 1], ] - grid$locs_grid[[1]][MRA$D$ind[, 2], ])
#     MRA$D$diffsx <- diffs$x
#     ddistx <- spam(MRA$D[c("ind", "diffsx")], nrow = MRA$D$da[1], ncol = MRA$D$da[2])
#     MRA$D$diffsy <- diffs$y
#     ddisty <- spam(MRA$D[c("ind", "diffsy")], nrow = MRA$D$da[1], ncol = MRA$D$da[2])
#     return(list(ddistx=ddistx, ddisty=ddisty))
# }

update_deep_mra <- function(y, locs, grid, MRA, MRA1, MRA2,
                            alpha,
                            alpha_x1, alpha_y1,
                            alpha_x2, alpha_y2,
                            Q,
                            Q1,
                            Q2,
                            i,
                            m, v,
                            learn_rate) {

    N <- length(y)
    # first layer
    delta <- drop(- 1/N * (y - MRA$W %*% alpha))

    # second layer
    delta_x1 <- delta * ((MRA$dW * MRA$ddistx) %*% alpha_x1)# check these
    delta_y1 <- delta * ((MRA$dW * MRA$ddisty) %*% alpha_y1)# check these

    # third layer
    delta_x2 <- (delta_x1 + delta_y1) * ((MRA1$dW * MRA1$ddistx) %*% alpha_x2)# check these
    delta_y2 <- (delta_x1 + delta_y1) * ((MRA1$dW * MRA1$ddisty) %*% alpha_y2)# check these

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
                    alpha = NULL,
                    alpha_x1 = NULL,
                    alpha_y1 = NULL,
                    alpha_x2 = NULL,
                    alpha_y2 = NULL,
                    learn_rate, n_iter=500,
                    n_message = 50) {


    N <- length(y)


    # initialize the MRA parameters here
    MRA2 <- eval_basis(locs, grid)
    Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=rep(0.9, length(MRA2$n_dims)))
    if (length(MRA2$n_dims) > 1) {
        Q2 <- do.call(bdiag.spam, Q2)
    }
    class(Q2) <- "spam"
    if (is.null(alpha_x2)) {
        alpha_x2 <- rnorm(ncol(MRA2$W), 0, 0.1)
    }
    if (is.null(alpha_y2)) {
        alpha_y2 <- rnorm(ncol(MRA2$W), 0, 0.1)
    }


    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid)
    Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=rep(0.9, length(MRA1$n_dims)))
    if (length(MRA1$n_dims) > 1) {
        Q1 <- do.call(bdiag.spam, Q1)
    }
    class(Q1) <- "spam"

    if (is.null(alpha_x1)) {
        alpha_x1 <- rnorm(ncol(MRA1$W), 0, 0.1)
    }
    if (is.null(alpha_y1)) {
        alpha_y1 <- rnorm(ncol(MRA1$W), 0, 0.1)
    }

    MRA <- eval_basis(cbind(MRA1$W %*% alpha_x1, MRA1$W %*% alpha_y1), grid)
    Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=rep(0.9, length(MRA1$n_dims)))
    if (length(MRA$n_dims) > 1) {
        Q <- do.call(bdiag.spam, Q)
    }
    class(Q) <- "spam"

    if (is.null(alpha)) {
        alpha <- rnorm(ncol(MRA$W), 0, 0.1)
    }

    # initialize the loss
    loss <- rep(NA, n_iter)

    # initialize adam optimization
    m <- vector(mode='list', length = 5) # alpha, alpha_x1, and alpha_y1
    v <- vector(mode='list', length = 5) # alpha, alpha_x1, and alpha_y1
    m[[1]] <- rep(0, length(alpha))
    m[[2]] <- rep(0, length(alpha_x1))
    m[[3]] <- rep(0, length(alpha_y1))
    m[[4]] <- rep(0, length(alpha_x2))
    m[[5]] <- rep(0, length(alpha_y2))
    v[[1]] <- rep(0, length(alpha))
    v[[2]] <- rep(0, length(alpha_x1))
    v[[3]] <- rep(0, length(alpha_y1))
    v[[4]] <- rep(0, length(alpha_x2))
    v[[5]] <- rep(0, length(alpha_y2))

    for (i in 1:n_iter) {


        pars <- update_deep_mra(y, locs, grid, MRA, MRA1, MRA2,
                                alpha,
                                alpha_x1, alpha_y1,
                                alpha_x2, alpha_y2,
                                # alpha_x2, alpha_y2,
                                Q, Q1, Q2,
                                i,
                                m, v,
                                learn_rate)
        alpha <- pars$alpha
        alpha_x1 <- pars$alpha_x1
        alpha_y1 <- pars$alpha_y1
        alpha_x2 <- pars$alpha_x2
        alpha_y2 <- pars$alpha_y2
        MRA <- pars$MRA
        MRA1 <- pars$MRA1
        m <- pars$m
        v <- pars$v

        loss[i] <- 1 / N * sum((y - MRA$W %*% alpha)^2)
        if (i %% n_message == 0) {
            message("iteration i = ", i, " loss = ", loss[i])
        }
    }

    return(list(alpha = alpha, alpha_x1 = alpha_x1, alpha_y1 = alpha_y1,
                alpha_x2 = alpha_x2, alpha_y2 = alpha_y2,
                MRA = MRA, MRA1 = MRA1, MRA2 = MRA2, loss = loss))
}



