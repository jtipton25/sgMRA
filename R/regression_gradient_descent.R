#' Fit a regression gradient descent
#'
#' @param y The observed data.
#' @param X The design matrix for the fixed effects.
#' @param W The basis matrix for the random effects.
#' @param inits If NULL, random initialization is used. Otherwise it must
#' be a numeric vector of length equal to the total columns of X and W.
#' @param threshold The threshold for determining convergence.
#' @param alpha The learning rate.
#' @param num_iters The number of iterations to run the gradient descent algorithm.
#' @param print_every The number of iterations between printing the output.
#' @param minibatch_size If NULL, the full data are used. If set as a positive numeric value,
#' this is the number of samples used in the minibatch.
#'
#' @return The output of the gradient descent algorithm
#' @export
#'
#' @import spam
#' @import Matrix
#'
regression_gradient_descent <- function(y, X, W,
                                        inits = NULL, threshold = 0.01, learn_rate = 0.001,
                                        num_iters = 10, print_every = 1,
                                        minibatch_size = NULL){


    message("Initializing parameters")
    U <- cbind(X, W)
    tU <- t(U)
    tUU <- NULL
    tUy <- NULL
    if (is.null(minibatch_size)) {
        tUU <- tU %*% U
        tUy <- tU %*% y
    }

    # initialize algorithm
    loss <- rep(0, num_iters)
    beta_save <- matrix(0, nrow = num_iters, ncol = ncol(U))
    colnames(beta_save) <- paste0("beta[", 0:(ncol(U)-1),"]")
    if (is.null(inits)) {
        inits <- rep(0, ncol(U))
    }
    beta_old <- inits

    # initializing adam
    m <- list(rep(0, length(beta_old)))
    v <- list(rep(1, length(beta_old)))

    message("Starting gradient descent")
    loss[1] <- target_fun(y, U, beta_old)
    message("Initial loss = ", loss[1])

    beta_save[1, ] <- beta_old
    # first iteration
    if (is.null(minibatch_size)) {
        grad <- gradient_fun(y, tUy, tUU, beta_old)
        adam_out <- adam(1, list(grad), m, v)
        beta_new <- beta_old - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
    } else {
        idx <- sample(1:length(y), minibatch_size)
        y_idx <- y[idx]
        tU_idx <- tU[, idx]
        tUy_idx <- tU_idx %*% y_idx
        # tUU_idx <- tU_idx %*% t(tU_idx)
        # faster to subset tU to get tU_idx then transpose to U_idx than to subset U directly
        grad <- gradient_fun_mini(y_idx, tUy_idx, t(tU_idx), tU_idx, beta_old)
        adam_out <- adam(1, list(grad), m, v)
        beta_new <- beta_old - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
    }

    i <- 2
    print_every <- ifelse(print_every > 0, print_every, num_iters)
    # gradient descent loop
    while( any(abs(beta_new - beta_old) > threshold) & i <= num_iters){


        beta_old <- beta_new
        if (is.null(minibatch_size)) {
            grad <- gradient_fun(y, tUy, tUU, beta_new)
            adam_out <- adam(i, list(grad), m, v)
            beta_new <- beta_new - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
        } else {
            idx <- sample(1:length(y), minibatch_size)
            y_idx <- y[idx]
            tU_idx <- tU[, idx]
            tUy_idx <- tU_idx %*% y_idx
            # tUU_idx <- tU_idx %*% t(tU_idx)
            # faster to subset tU to get tU_idx then transpose to U_idx than to subset U directly
            grad <- gradient_fun_mini(y_idx, tUy_idx, t(tU_idx), tU_idx, beta_new)
            adam_out <- adam(i, list(grad), m, v)
            beta_new <- beta_new - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
        }
        loss[i] <- target_fun(y, U, beta_new)
        # save the estimates on the original data scale (not normalized scale)
        # beta_save[i, ] <- c(beta_new[1] * y_sd + y_mean - sum(beta_new[-1][1:(ncol(X)-1)] * X_mean / X_sd) * y_sd,
        #                     beta_new[-1] * y_sd / X_sd)
        beta_save[i, ] <- as.numeric(beta_new)
        if(i %% print_every == 0) {
            message("iteration = ", i, ", loss = ", loss[i])
        }
        i <- i+1
    }

    cbind(as.data.frame(beta_save[1:(i-1), ]),
          loss = loss[1:(i-1)],
          iteration = 1:(i-1))
}
