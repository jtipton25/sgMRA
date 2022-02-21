deep_mra_mh <- function(y, locs, n_mcmc = 1000,
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

  source(here::here("R", "update-tuning.R"))
  N <- length(y)

  grid <- make_grid(locs, M = M, n_coarse_grid = n_coarse_grid)

  # initialize the algorithm ----
  # construct the first layer
  # W1 <- mra_wendland_2d(as.matrix(locs), M=M, n_coarse_grid = n_coarse_grid)$W
  W1 <- eval_basis(as.matrix(locs), grid)$W
  Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=0.9)
  class(Q1) <- "spam"
  if (is.null(alpha_x1)) {
    alpha_x1 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
  }
  if (is.null(alpha_y1)) {
    alpha_y1 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
  }

  # # construct the second layer
  # W2 <- mra_wendland_2d(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), M=M, n_coarse_grid = n_coarse_grid)$W
  W2 <- eval_basis(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), grid)$W
  Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=0.9)
  class(Q2) <- "spam"
  if (is.null(alpha_x2)) {
    alpha_x2 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
  }
  if (is.null(alpha_y2)) {
    alpha_y2 <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q2)), Q = Q2))
  }

  # construct the final layer
  # W <- mra_wendland_2d(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), M=M, n_coarse_grid = n_coarse_grid)$W
  W <- eval_basis(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), grid)$W
  Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=0.9)
  class(Q) <- "spam"
  if (is.null(alpha)) {
    alpha <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q)), Q = Q))
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

  alpha_x1_tune <- rep(0.01, length(alpha_x1))
  alpha_x1_accept <- rep(0, length(alpha_x1))
  alpha_x1_accept_batch <- rep(0, length(alpha_x1))

  alpha_y1_tune <- rep(0.01, length(alpha_y1))
  alpha_y1_accept <- rep(0, length(alpha_y1))
  alpha_y1_accept_batch <- rep(0, length(alpha_y1))

  alpha_x2_tune <- rep(0.01, length(alpha_x2))
  alpha_x2_accept <- rep(0, length(alpha_x2))
  alpha_x2_accept_batch <- rep(0, length(alpha_x2))

  alpha_y2_tune <- rep(0.01, length(alpha_y2))
  alpha_y2_accept <- rep(0, length(alpha_y2))
  alpha_y2_accept_batch <- rep(0, length(alpha_y2))

  alpha_tune    <- rep(0.01, length(alpha))
  alpha_accept <- rep(0, length(alpha))
  alpha_accept_batch <- rep(0, length(alpha))

  # setup the ess pars ----
  # pars <- list(
  #     y             = y,
  #     alpha_x1      = alpha_x1,
  #     alpha_y1      = alpha_y1,
  #     W1            = W1,
  #     alpha_x2      = alpha_x2,
  #     alpha_y2      = alpha_y2,
  #     W2            = W2,
  #     alpha         = alpha,
  #     W             = W,
  #     sigma         = sigma,
  #     M             = M,
  #     n_coarse_grid = n_coarse_grid,
  #     num_calls_ess = 0,
  #     verbose = TRUE)


  message("Starting MCMC, will run for ", n_mcmc, " iterations")
  for (k in 1:n_mcmc) {
    if (k %% n_message == 0) {
      message("On iteration ", k, " out of ", n_mcmc)
    }

    # update alpha_x1 ----
    if (sample_alpha_x1) {
      message("sampling alpha_x1")
      for (j in 1:length(alpha_x1)) {
        alpha_x1_prop <- alpha_x1
        alpha_x1_prop[j] <- rnorm(1, alpha_x1[j], alpha_x1_tune[j])
        # W2_prop <- mra_wendland_2d(cbind(W1 %*% alpha_x1_prop, W1 %*% alpha_y1), M=M, n_coarse_grid = n_coarse_grid)$W
        W2_prop <- eval_basis(cbind(W1 %*% alpha_x1_prop, W1 %*% alpha_y1), grid)$W
        # W_prop  <- mra_wendland_2d(cbind(W2_prop %*% alpha_x2, W2_prop %*% alpha_y2), M=M, n_coarse_grid = n_coarse_grid)$W
        W_prop  <- eval_basis(cbind(W2_prop %*% alpha_x2, W2_prop %*% alpha_y2), grid)$W
        mh1 <- sum(dnorm(y, W_prop %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_x1_prop * Q1 %*% alpha_x1_prop)
        mh2 <- sum(dnorm(y, W %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_x1 * Q1 %*% alpha_x1)
        mh <- exp(mh1 - mh2)
        if (mh > runif(1)) {
          alpha_x1 <- alpha_x1_prop
          W2 <- W2_prop
          W  <- W_prop
          alpha_x1_accept_batch <- alpha_x1_accept_batch + 1 / 50
          alpha_x1_accept <- alpha_x1_accept + 1 / n_mcmc
        }
      }

      if (k %% 50 == 0) {
        out_tune <- update_tuning_vec(k, alpha_x1_accept_batch, alpha_x1_tune)
        alpha_x1_accept_batch <- out_tune$accept
        alpha_x1_tune <- out_tune$tune
      }
    }

    # update alpha_y1 ----
    if (sample_alpha_y1) {

      message("sampling alpha_y1")
      for (j in 1:length(alpha_y1)) {
        alpha_y1_prop <- alpha_y1
        alpha_y1_prop[j] <- rnorm(1, alpha_y1[j], alpha_y1_tune[j])
        # W2_prop <- mra_wendland_2d(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1_prop), M=M, n_coarse_grid = n_coarse_grid)$W
        W2_prop <- eval_basis(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1_prop), grid)$W
        # W_prop  <- mra_wendland_2d(cbind(W2_prop %*% alpha_x2, W2_prop %*% alpha_y2), M=M, n_coarse_grid = n_coarse_grid)$W
        W_prop  <- eval_basis(cbind(W2_prop %*% alpha_x2, W2_prop %*% alpha_y2), grid)$W
        mh1 <- sum(dnorm(y, W_prop %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_y1_prop * Q1 %*% alpha_y1_prop)
        mh2 <- sum(dnorm(y, W %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_y1 * Q1 %*% alpha_y1)
        mh <- exp(mh1 - mh2)
        if (mh > runif(1)) {
          alpha_y1 <- alpha_y1_prop
          W2 <- W2_prop
          W  <- W_prop
          alpha_y1_accept_batch <- alpha_y1_accept_batch + 1 / 50
          alpha_y1_accept <- alpha_y1_accept + 1 / n_mcmc
        }
      }

      if (k %% 50 == 0) {
        out_tune <- update_tuning_vec(k, alpha_y1_accept_batch, alpha_y1_tune)
        alpha_y1_accept_batch <- out_tune$accept
        alpha_y1_tune <- out_tune$tune
      }
    }

    # update alpha_x2 ----
    if (sample_alpha_x2) {
      message("sampling alpha_x2")
      for (j in 1:length(alpha_x2)) {
        alpha_x2_prop <- alpha_x2
        alpha_x2_prop[j] <- rnorm(1, alpha_x2[j], alpha_x2_tune[j])
        # W_prop  <- mra_wendland_2d(cbind(W2 %*% alpha_x2_prop, W2 %*% alpha_y2), M=M, n_coarse_grid = n_coarse_grid)$W
        W_prop  <- eval_basis(cbind(W2 %*% alpha_x2_prop, W2 %*% alpha_y2), grid)$W
        mh1 <- sum(dnorm(y, W_prop %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_x2_prop * Q2 %*% alpha_x2_prop)
        mh2 <- sum(dnorm(y, W %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_x2 * Q2 %*% alpha_x2)
        mh <- exp(mh1 - mh2)
        if (mh > runif(1)) {
          alpha_x2 <- alpha_x2_prop
          W  <- W_prop
          alpha_x2_accept_batch <- alpha_x2_accept_batch + 1 / 50
          alpha_x2_accept <- alpha_x2_accept + 1 / n_mcmc
        }
      }

      if (k %% 50 == 0) {
        out_tune <- update_tuning_vec(k, alpha_x2_accept_batch, alpha_x2_tune)
        alpha_x2_accept_batch <- out_tune$accept
        alpha_x2_tune <- out_tune$tune
      }
    }

    # update alpha_y2 ----
    if (sample_alpha_y2) {
      message("sampling alpha_y2")
      for (j in 1:length(alpha_y2)) {
        alpha_y2_prop <- alpha_y2
        alpha_y2_prop[j] <- rnorm(1, alpha_y2[j], alpha_y2_tune[j])
        # W_prop  <- mra_wendland_2d(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2_prop), M=M, n_coarse_grid = n_coarse_grid)$W
        W_prop  <- eval_basis(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2_prop), grid)$W
        mh1 <- sum(dnorm(y, W_prop %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_y2_prop * Q2 %*% alpha_y2_prop)
        mh2 <- sum(dnorm(y, W %*% alpha, sigma, log=TRUE)) -
          0.5 * sum(alpha_y2 * Q2 %*% alpha_y2)
        mh <- exp(mh1 - mh2)
        if (mh > runif(1)) {
          alpha_y2 <- alpha_y2_prop
          W  <- W_prop
          alpha_y2_accept_batch <- alpha_y2_accept_batch + 1 / 50
          alpha_y2_accept <- alpha_y2_accept + 1 / n_mcmc
        }
      }

      if (k %% 50 == 0) {
        out_tune <- update_tuning_vec(k, alpha_y2_accept_batch, alpha_y2_tune)
        alpha_y2_accept_batch <- out_tune$accept
        alpha_y2_tune <- out_tune$tune
      }
    }

    # update alpha ----
    if (sample_alpha) {
      message("sampling alpha")
      # sample using MH

      # alpha_prop <- rnorm(length(alpha), alpha, alpha_tune)
      # mh1 <- sum(dnorm(y, W %*% alpha_prop, sigma, log=TRUE)) -
      #     0.5 * sum(alpha_prop * Q %*% alpha_prop)
      # mh2 <- sum(dnorm(y, W %*% alpha, sigma, log=TRUE)) -
      #     0.5 * sum(alpha * Q %*% alpha)
      # mh <- exp(mh1 - mh2)
      # if (mh > runif(1)) {
      #     alpha <- alpha_prop
      #     alpha_accept_batch <- alpha_accept_batch + 1 / 50
      #     alpha_accept <- alpha_accept + 1 / n_mcmc
      # }
      #
      # if (k %% 50 == 0) {
      #     out_tune <- update_tuning(k, alpha_accept_batch, alpha_tune)
      #     alpha_accept_batch <- out_tune$accept
      #     alpha_tune <- out_tune$tune
      # }

      # sample using Gibbs update
      A_alpha <- 1 / sigma^2 * t(W) %*% W + Q
      b_alpha <- 1 / sigma^2 * t(W) %*% y
      alpha   <- as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha))

    }

    # update sigma ----
    message("sampling sigma")
    sigma <- sqrt(1 / rgamma(1, 1 + 0.5 * N, 1 + 0.5 * sum((y - W %*% alpha)^2)))

    message("sigma = ", sigma)

    # save the parameters
    alpha_x1_save[k, ] <- alpha_x1
    alpha_y1_save[k, ] <- alpha_y1
    # W1_save[[k]]       <- W1
    alpha_x2_save[k, ] <- alpha_x2
    alpha_y2_save[k, ] <- alpha_y2
    # W2_save[[k]]       <- W2
    alpha_save[k, ]    <- alpha
    # W_save[[k]]        <- W
    Walpha_save[k, ]   <- W %*% alpha
    sigma_save[k]      <- sigma

  }

  message("Acceptance rate for alpha_x1 = ", alpha_x1_accept)
  message("Acceptance rate for alpha_y1 = ", alpha_y1_accept)
  message("Acceptance rate for alpha_x2 = ", alpha_x2_accept)
  message("Acceptance rate for alpha_y2 = ", alpha_y2_accept)
  message("Acceptance rate for alpha = ", alpha_accept)

  return(list(alpha_x1 = alpha_x1_save,
              alpha_y1 = alpha_y1_save,
              # W1       = W1_save,
              alpha_x2 = alpha_x2_save,
              alpha_y2 = alpha_y2_save,
              # W2       = W2_save,
              alpha    = alpha_save,
              # W        = W_save,
              Walpha   = Walpha_save,
              sigma    = sigma_save))
}
