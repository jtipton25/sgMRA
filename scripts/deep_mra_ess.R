deep_mra_ess <- function(y, locs, n_mcmc = 1000, M = 1, n_coarse_grid = 20,
                         n_message = 50) {

    library(spam)
    library(Matrix)
    library(igraph)
    library(tidyverse)
    library(BayesMRA)

    N <- length(y)

    # helper functions for ess
    update_pars_alpha_x1 <- function(alpha_x1_prop, pars) {
        pars$alpha_x1 <- alpha_x1_prop
        pars$W2 <- mra_wendland_2d(cbind(pars$W1 %*% pars$alpha_x1, pars$W1 %*% pars$alpha_y1), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        pars$W  <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        return(pars)
    }

    update_pars_alpha_y1 <- function(alpha_y1_prop, pars) {
        pars$alpha_y1 <- alpha_y1_prop
        pars$W2 <- mra_wendland_2d(cbind(pars$W1 %*% pars$alpha_x1, pars$W1 %*% pars$alpha_y1), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        pars$W  <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        return(pars)
    }

    update_pars_alpha_x2 <- function(alpha_x2_prop, pars) {
        pars$alpha_x2 <- alpha_x2_prop
        pars$W <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        return(pars)
    }

    update_pars_alpha_y2 <- function(alpha_y2_prop, pars) {
        pars$alpha_y2 <- alpha_y2_prop
        pars$W <- mra_wendland_2d(cbind(pars$W2 %*% pars$alpha_x2, pars$W2 %*% pars$alpha_y2), M=pars$M, n_coarse_grid = pars$n_coarse_grid)$W
        return(pars)
    }

    update_pars_alpha <- function(alpha_prop, pars) {
        pars$alpha <- alpha_prop
        return(pars)
    }

    log_like <- function(pars) {
        sum(dnorm(pars$y, pars$W %*% pars$alpha, pars$sigma, log=TRUE))
    }

    # initialize the algorithm ----
    # construct the first layer
    W1 <- mra_wendland_2d(as.matrix(locs), M=M, n_coarse_grid = n_coarse_grid)$W
    Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=0.9)
    class(Q1) <- "spam"
    alpha_x1 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
    alpha_y1 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
    # construct the second layer
    W2 <- mra_wendland_2d(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), M=M, n_coarse_grid = n_coarse_grid)$W
    Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=0.9)
    class(Q2) <- "spam"
    alpha_x2 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
    alpha_y2 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))

    # construct the final layer
    W <- mra_wendland_2d(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), M=M, n_coarse_grid = n_coarse_grid)$W
    Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=0.9)
    class(Q) <- "spam"
    alpha <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q)), Q = Q))

    sigma <- max(min(rgamma(1, 1, 1), 5), 0.05)

    alpha_x1_save <- matrix(0, n_mcmc, nrow(Q1))
    alpha_y1_save <- matrix(0, n_mcmc, nrow(Q1))
    W1_save       <- vector('list', n_mcmc)
    alpha_x2_save <- matrix(0, n_mcmc, nrow(Q2))
    alpha_y2_save <- matrix(0, n_mcmc, nrow(Q2))
    W2_save       <- vector('list', n_mcmc)
    alpha_save    <- matrix(0, n_mcmc, nrow(Q))
    W_save        <- vector('list', n_mcmc)
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
        message("sampling alpha_x1")
        alpha_x1_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
        pars <- ess(current=alpha_x1, prior=alpha_x1_prop,
                    prior_mean=rep(0, length(alpha_x1)), pars=pars,
                    log_like_fun=log_like, update_pars_fun=update_pars_alpha_x1)

        # update alpha_y1 ----
        message("sampling alpha_y1")
        alpha_y1_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
        pars <- ess(current=alpha_y1, prior=alpha_y1_prop,
                    prior_mean=rep(0, length(alpha_y1)), pars=pars,
                    log_like_fun=log_like, update_pars_fun=update_pars_alpha_y1)

        # update alpha_x2 ----
        message("sampling alpha_x2")
        alpha_x2_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
        pars <- ess(current=alpha_x2, prior=alpha_x2_prop,
                    prior_mean=rep(0, length(alpha_x2)), pars=pars,
                    log_like_fun=log_like, update_pars_fun=update_pars_alpha_x2)

        # update alpha_y2 ----
        message("sampling alpha_y2")
        alpha_y2_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
        pars <- ess(current=alpha_y2, prior=alpha_y2_prop,
                    prior_mean=rep(0, length(alpha_y2)), pars=pars,
                    log_like_fun=log_like, update_pars_fun=update_pars_alpha_y2)

        # update alpha ----
        message("sampling alpha")
        alpha_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q)), Q = Q))
        pars <- ess(current=alpha, prior=alpha_prop,
                    prior_mean=rep(0, length(alpha)), pars=pars,
                    log_like_fun=log_like, update_pars_fun=update_pars_alpha)

        # update sigma ----
        message("sampling sigma")
        pars$sigma <- sqrt(1 / rgamma(1, 1 + 0.5 * N, 1 + 0.5 * sum((pars$y - pars$W %*% pars$alpha)^2)))

        message("sigma = ", pars$sigma)

        # save the parameters
        alpha_x1_save[k, ] <- pars$alpha_x1
        alpha_y1_save[k, ] <- pars$alpha_y1
        W1_save[[k]]       <- pars$W1
        alpha_x2_save[k, ] <- pars$alpha_x2
        alpha_y2_save[k, ] <- pars$alpha_y2
        W2_save[[k]]       <- pars$W2
        alpha_save[k, ]    <- pars$alpha
        W_save[[k]]        <- pars$W
        sigma_save[k  ]    <- pars$sigma

    }

    return(list(alpha_x1 = alpha_x1_save,
                alpha_y1 = alpha_y1_save,
                W1       = W1_save,
                alpha_x2 = alpha_x2_save,
                alpha_y2 = alpha_y2_save,
                W2       = W2_save,
                alpha    = alpha_save,
                W        = W_save,
                sigma    = sigma$save))
}
