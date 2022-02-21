deep_mra_ess <- function(y, locs, n_mcmc = 1000,
                         alpha_x1 = NULL,
                         sample_alpha_x1 = TRUE,
                         alpha_y1 = NULL,
                         sample_alpha_y1 = TRUE,
                         alpha_x2 = NULL,
                         sample_alpha_x2 = TRUE,
                         alpha_y2 = NULL,
                         sample_alpha_y2 = TRUE,
                         alpha = NULL,
                         sample_alpha = TRUE,
                         M = 1, n_coarse_grid = 20,
                         n_message = 50) {

    library(spam)
    library(Matrix)
    library(igraph)
    library(tidyverse)
    library(BayesMRA)

    N <- length(y)

    grid <- make_grid(locs, M = M, n_coarse_grid = n_coarse_grid)

    # helper functions for ess
    update_pars_alpha_x1 <- function(alpha_x1_prop, pars) {
        pars$alpha_x1 <- alpha_x1_prop
        # pars$W2 <- mra_wendland_2d(cbind(pars$W1 %*% pars$alpha_x1, pars$W1 %*% pars$alpha_y1), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        # pars$W  <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        pars$W2 <- eval_basis(cbind(pars$W1 %*% pars$alpha_x1, pars$W1 %*% pars$alpha_y1), pars$grid)$W
        pars$W  <- eval_basis(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), pars$grid)$W
        return(pars)
    }

    update_pars_alpha_y1 <- function(alpha_y1_prop, pars) {
        pars$alpha_y1 <- alpha_y1_prop
        # pars$W2 <- mra_wendland_2d(cbind(pars$W1 %*% pars$alpha_x1, pars$W1 %*% pars$alpha_y1), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        # pars$W  <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        pars$W2 <- eval_basis(cbind(pars$W1 %*% pars$alpha_x1, pars$W1 %*% pars$alpha_y1), pars$grid)$W
        pars$W  <- eval_basis(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), pars$grid)$W
        return(pars)
    }

    update_pars_alpha_x2 <- function(alpha_x2_prop, pars) {
        pars$alpha_x2 <- alpha_x2_prop
        # pars$W <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        pars$W  <- eval_basis(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), pars$grid)$W
        return(pars)
    }

    update_pars_alpha_y2 <- function(alpha_y2_prop, pars) {
        pars$alpha_y2 <- alpha_y2_prop
        # pars$W <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        pars$W  <- eval_basis(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), pars$grid)$W
        return(pars)
    }

    update_pars_alpha <- function(alpha_prop, pars) {
        pars[['alpha']] <- alpha_prop
        return(pars)
    }

    log_like <- function(pars) {
        return(sum(dnorm(pars$y, pars$W %*% pars$alpha, pars$sigma, log=TRUE)))
    }

    # initialize the algorithm ----
    # construct the first layer
    # W1 <- mra_wendland_2d(as.matrix(locs), M=M, n_coarse_grid = n_coarse_grid)$
    W1 <- eval_basis(as.matrix(locs), grid)$W
    Q1 <- make_Q_alpha_2d(sqrt(ncol(W1)), phi=0.9)
    class(Q1) <- "spam"

    if (is.null(alpha_x1)) {
        # alpha_x1 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
        alpha_x1 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=rep(1, nrow(W1)) %*% W1, a=0))
    }
    if (is.null(alpha_y1)) {
        # alpha_y1 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
        alpha_y1 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=rep(1, nrow(W1)) %*% W1, a=0))
    }

    # construct the second layer
    # W2 <- mra_wendland_2d(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), M=M, n_coarse_grid = n_coarse_grid)$W
    W2 <- eval_basis(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), grid)$W
    Q2 <- make_Q_alpha_2d(sqrt(ncol(W2)), phi=0.9)
    class(Q2) <- "spam"

    if (is.null(alpha_x2)) {
        # alpha_x2 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
        alpha_x2 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=rep(1, nrow(W2)) %*% W2, a=0))
    }
    if (is.null(alpha_y2)) {
        # alpha_y2 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
        alpha_y2 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=rep(1, nrow(W2)) %*% W2, a=0))
    }
    # construct the final layer
    # W <- mra_wendland_2d(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), M=M, n_coarse_grid = n_coarse_grid)$W
    W <- eval_basis(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), grid)$W
    Q <- make_Q_alpha_2d(sqrt(ncol(W)), phi=0.9)
    class(Q) <- "spam"
    if (is.null(alpha)) {
        # alpha <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q)), Q = Q))
        alpha <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q)), Q = Q, A=rep(1, nrow(W)) %*% W, a=0))
    }

    sigma <- max(min(rgamma(1, 1, 1), 5), 0.05)

    alpha_x1_save <- matrix(0, n_mcmc, nrow(Q1))
    alpha_y1_save <- matrix(0, n_mcmc, nrow(Q1))
    # W1_save       <- vector('list', n_mcmc)
    alpha_x2_save <- matrix(0, n_mcmc, nrow(Q2))
    alpha_y2_save <- matrix(0, n_mcmc, nrow(Q2))
    # W2_save       <- vector('list', n_mcmc)
    alpha_save    <- matrix(0, n_mcmc, nrow(Q))
    # W_save        <- vector('list', n_mcmc)
    # Walpha1_save  <- matrix(0, n_mcmc, N)
    # Walpha2_save  <- matrix(0, n_mcmc, N)
    Walpha_save   <- matrix(0, n_mcmc, N)
    sigma_save    <- rep(0, n_mcmc)

    # setup the ess pars ----
    pars <- list(
        y             = y,
        alpha_x1      = alpha_x1,
        alpha_y1      = alpha_y1,
        W1            = W1,
        alpha_x2      = alpha_x2,
        alpha_y2      = alpha_y2,
        W2            = W2,
        alpha         = alpha,
        W             = W,
        sigma         = sigma,
        grid          = grid,
        M             = M,
        n_coarse_grid = n_coarse_grid,
        num_calls_ess = 0,
        verbose = TRUE)


    message("Starting MCMC, will run for ", n_mcmc, " iterations")
    for (k in 1:n_mcmc) {
        if (k %% n_message == 0) {
            message("On iteration ", k, " out of ", n_mcmc)
        }

        # update alpha_x1 ----
        if (sample_alpha_x1) {
            message("sampling alpha_x1")
            alpha_x1_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
            # alpha_x1_prop <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=rep(1, nrow(W1)) %*% W1, a=0))
            pars <- ess(current=pars$alpha_x1, prior=alpha_x1_prop,
                        prior_mean=rep(0, length(pars$alpha_x1)), pars=pars,
                        log_like_fun=log_like, update_pars_fun=update_pars_alpha_x1)
        }

        # update alpha_y1 ----
        if (sample_alpha_y1) {
            message("sampling alpha_y1")
            alpha_y1_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
            # alpha_y1_prop <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=rep(1, nrow(W1)) %*% W1, a=0))
            pars <- ess(current=pars$alpha_y1, prior=alpha_y1_prop,
                        prior_mean=rep(0, length(pars$alpha_y1)), pars=pars,
                        log_like_fun=log_like, update_pars_fun=update_pars_alpha_y1)
        }

        # update alpha_x2 ----
        if (sample_alpha_x2) {
            message("sampling alpha_x2")
            alpha_x2_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
            # alpha_x2_prop <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=rep(1, nrow(W2)) %*% W2, a=0))
            # alpha_x2_prop <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=t(rep(1, ncol(W2))), a=0))
            pars <- ess(current=pars$alpha_x2, prior=alpha_x2_prop,
                        prior_mean=rep(0, length(pars$alpha_x2)), pars=pars,
                        log_like_fun=log_like, update_pars_fun=update_pars_alpha_x2)
        }

        # update alpha_y2 ----
        if (sample_alpha_y2) {
            message("sampling alpha_y2")
            alpha_y2_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
            # alpha_y2_prop <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=rep(1, nrow(W2)) %*% W2, a=0))
            # alpha_y2_prop <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=t(rep(1, ncol(W2))), a=0))
            pars <- ess(current=pars$alpha_y2, prior=alpha_y2_prop,
                        prior_mean=rep(0, length(pars$alpha_y2)), pars=pars,
                        log_like_fun=log_like, update_pars_fun=update_pars_alpha_y2)
        }

        # update alpha ----
        if (sample_alpha) {
            message("sampling alpha")
            # alpha_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q)), Q = Q, A=rep(1, nrow(pars$W)) %*% pars$W, a=0))
            # pars <- ess(current=pars$alpha, prior=alpha_prop,
            #             prior_mean=rep(0, length(pars$alpha)), pars=pars,
            #             log_like_fun=log_like, update_pars_fun=update_pars_alpha)
            A_alpha <- 1 / pars$sigma^2 * t(pars$W) %*% pars$W + Q
            b_alpha <- 1 / pars$sigma^2 * t(pars$W) %*% pars$y
            # alpha   <- as.vector(rmvnorm.canonical(1, b_alpha, A_alpha))
            alpha   <- as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, A=rep(1, nrow(pars$W)) %*% pars$W, a=0))
            # alpha   <- as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, A=t(rep(1, ncol(W))), a=0))
            pars$alpha <- alpha
        }

        # update sigma ----
        message("sampling sigma")
        pars$sigma <- sqrt(1 / rgamma(1, 1 + 0.5 * N, 1 + 0.5 * sum((pars$y - pars$W %*% pars$alpha)^2)))

        message("sigma = ", pars$sigma)

        # save the parameters
        alpha_x1_save[k, ] <- pars$alpha_x1
        alpha_y1_save[k, ] <- pars$alpha_y1
        # W1_save[[k]]       <- pars$W1
        alpha_x2_save[k, ] <- pars$alpha_x2
        alpha_y2_save[k, ] <- pars$alpha_y2
        # W2_save[[k]]       <- pars$W2
        alpha_save[k, ]    <- pars$alpha
        # W_save[[k]]        <- pars$W
        Walpha_save[k, ]   <- pars$W %*% pars$alpha
        sigma_save[k]      <- pars$sigma

    }

    return(list(alpha_x1 = alpha_x1_save,
                alpha_y1 = alpha_y1_save,
                # W1       = W1_save,
                alpha_x2 = alpha_x2_save,
                alpha_y2 = alpha_y2_save,
                # W2       = W2_save,
                alpha    = alpha_save,
                # W       = W_save,
                # Walpha1  = Walpha_save,
                # Walpha2  = Walpha_save,
                Walpha   = Walpha_save,
                grid     = grid,
                sigma    = sigma_save))
}
