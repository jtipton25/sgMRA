nonstationary_MRA_ess <- function(y,
                                  locs,
                                  grid,
                                  grid_kernels,
                                  radius = 0.5,
                                  trunc = 0.01,
                                  n_mcmc = 1000,
                                  alpha = NULL,
                                  alpha_a = NULL,
                                  alpha_b = NULL,
                                  sample_alpha = TRUE,
                                  sample_alpha_a = TRUE,
                                  sample_alpha_b = TRUE,
                                  n_message = 50) {

    # define helper functions

    update_pars_alpha_a <- function(alpha_a_prop, pars) {
        pars$alpha_a <- alpha_a_prop
        pars$a       <- pmax(drop(pars$W_kernels %*% pars$alpha_a), pars$trunc)
        for (i in 1:pars$N) {
            pars$D[i, ] <- sqrt(rowSums(pars$sq_devs[i, , ] * cbind(rep(pars$a[i], pars$N_grid), rep(pars$b[i], pars$N_grid))))
        }

        pars$W  <- as.spam(BayesMRA::wendland_basis(pars$D, pars$radius))
        return(pars)
    }



    update_pars_alpha_b <- function(alpha_b_prop, pars) {
        pars$alpha_b <- alpha_b_prop
        pars$b       <- pmax(drop(pars$W_kernels %*% pars$alpha_b), pars$trunc)
        for (i in 1:pars$N) {
            pars$D[i, ] <- sqrt(rowSums(pars$sq_devs[i, , ] * cbind(rep(pars$a[i], pars$N_grid), rep(pars$b[i], pars$N_grid))))
        }

        pars$W  <- as.spam(BayesMRA::wendland_basis(pars$D, pars$radius))
        return(pars)
    }


    log_like <- function(pars) {
        return(sum(dnorm(pars$y, pars$W %*% pars$alpha, pars$sigma, log=TRUE)))
    }


    N <- length(y)
    N_grid <- nrow(grid)
    N_grid_kernels <- nrow(grid_kernels)

    # initialize the kernel process
    D_kernels <- matrix(0, N, N_grid_kernels)
    for (i in 1:N) {
        D_kernels[i, ] <- sqrt(rowSums(sweep(grid_kernels, 2, FUN='-', locs[i, ])^2))
    }

    W_kernels <- as.spam(BayesMRA::wendland_basis(D_kernels, radius))
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


    D <- matrix(0, N, N_grid)
    for (i in 1:N) {
        D[i, ] <- sqrt(rowSums(sq_devs[i, , ] * cbind(rep(a[i], N_grid), rep(b[i], N_grid))))
    }

    W <- as.spam(BayesMRA::wendland_basis(D, radius))

    if (is.null(alpha)) {
        # add in penalty parameter later
        alpha <- rnorm(N_grid, 0, 0.1)
    }

    # setup save variables ----


    sigma <- max(min(rgamma(1, 1, 1), 5), 0.05)

    alpha_a_save <- matrix(0, n_mcmc, N_grid_kernels)
    alpha_b_save <- matrix(0, n_mcmc, N_grid_kernels)
    alpha_save    <- matrix(0, n_mcmc, nrow(Q))
    Walpha_save   <- matrix(0, n_mcmc, N)
    sigma_save    <- rep(0, n_mcmc)

    # setup the ess pars ----
    pars <- list(
        y             = y,
        alpha_a       = alpha_a,
        alpha_b       = alpha_b,
        a             = a,
        b             = b,
        alpha         = alpha,
        N             = N,
        N_grid        = N_grid,
        W             = W,
        W_kernels     = W_kernels,
        sq_devs       = sq_devs,
        D             = D,
        radius        = radius,
        trunc         = trunc,
        sigma         = sigma,
        num_calls_ess = 0,
        verbose = TRUE)


    message("Starting MCMC, will run for ", n_mcmc, " iterations")
    for (k in 1:n_mcmc) {
        if (k %% n_message == 0) {
            message("On iteration ", k, " out of ", n_mcmc)
        }

        # update alpha_a ----
        if (sample_alpha_a) {
            message("sampling alpha_a")
            alpha_a_prop <- rnorm(N_grid_kernels, 0, 0.1) #drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
            pars <- ess(current=pars$alpha_a, prior=alpha_a_prop,
                        prior_mean=rep(0, length(pars$alpha_a)), pars=pars,
                        log_like_fun=log_like, update_pars_fun=update_pars_alpha_a)
        }

        # update alpha_b ----
        if (sample_alpha_b) {
            message("sampling alpha_b")
            alpha_b_prop <- rnorm(N_grid_kernels, 0, 0.1) #drop(rmvnorm.prec(1, mu = rep(0, nrow(Q1)), Q = Q1))
            pars <- ess(current=pars$alpha_b, prior=alpha_b_prop,
                        prior_mean=rep(0, length(pars$alpha_b)), pars=pars,
                        log_like_fun=log_like, update_pars_fun=update_pars_alpha_b)
        }

        # update alpha ----
        if (sample_alpha) {
            message("sampling alpha")
            # alpha_prop <- drop(rmvnorm.prec(1, mu = rep(0, nrow(Q)), Q = Q, A=rep(1, nrow(pars$W)) %*% pars$W, a=0))
            # pars <- ess(current=pars$alpha, prior=alpha_prop,
            #             prior_mean=rep(0, length(pars$alpha)), pars=pars,
            #             log_like_fun=log_like, update_pars_fun=update_pars_alpha)
            A_alpha <- 1 / pars$sigma^2 * t(pars$W) %*% pars$W + diag(N_grid) / 100 #Q
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
        alpha_a_save[k, ] <- pars$alpha_a
        alpha_b_save[k, ] <- pars$alpha_b
        alpha_save[k, ]    <- pars$alpha
        Walpha_save[k, ]   <- pars$W %*% pars$alpha
        sigma_save[k]      <- pars$sigma

    }

    return(list(alpha_a = alpha_a_save,
                alpha_b = alpha_b_save,
                alpha   = alpha_save,
                Walpha  = Walpha_save,
                sigma   = sigma_save))
}
