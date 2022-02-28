# Deep BayesMRA
library(spam)
library(Matrix)
library(igraph)
library(tidyverse)
library(BayesMRA)
library(patchwork)

source("~/sgMRA/R/eval_basis.R")
Rcpp::sourceCpp("~/sgMRA/src/dist_near_cpp.cpp")
source("~/sgMRA/R/dwendland_basis.R")
source("~/sgMRA/R/adam.R")

penalized = TRUE
learn_rate=0.1

N <- 2^12

M <- 3
n_coarse_grid <- 30
# N <- 2^12
# M <- 1
# n_coarse_grid <- 80
source("~/sgMRA/R/sim-deep-mra.R")
dat_sim <- sim_deep_mra(N, M, n_coarse_grid, n_layers = 3, sigma = 0.1)


y=dat_sim$y
locs=dat_sim$locs
grid=dat_sim$grid
MRA=dat_sim$MRA[[1]]
MRA1=dat_sim$MRA[[2]]
MRA2=dat_sim$MRA[[3]]
W = MRA$W
W1 = MRA1$W
W2 = MRA2$W
dW = MRA$dW
dW1 = MRA1$dW
ddistx = MRA$ddistx
ddistx2 = MRA1$ddistx
ddisty = MRA$ddisty
ddisty2 = MRA1$ddisty

alpha=dat_sim$alpha
alpha_x1=dat_sim$alpha_x[[1]]
alpha_y1=dat_sim$alpha_y[[1]]
alpha_x2=dat_sim$alpha_x[[2]]
alpha_y2=dat_sim$alpha_y[[2]]

use_spam = TRUE


if (use_spam) {
    Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=rep(0.9, length(MRA2$n_dims)))
    if (length(MRA2$n_dims) > 1) {
        for (m in 1:M){
            Q2[[m]] <- 2^(2*(m-1)) * Q2[[m]]
        }
        Q2 <- do.call(bdiag.spam, Q2)
    }
    class(Q2) <- "spam"
    Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=rep(0.9, length(MRA1$n_dims)))
    if (length(MRA1$n_dims) > 1) {
        for (m in 1:M){
            Q1[[m]] <- 2^(2*(m-1)) * Q1[[m]]
        }
        Q1 <- do.call(bdiag.spam, Q1)
    }
    class(Q1) <- "spam"


    Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=rep(0.9, length(MRA1$n_dims)))
    if (length(MRA$n_dims) > 1) {
        for (m in 1:M){
            Q[[m]] <- 2^(2*(m-1)) * Q[[m]]
        }
        Q <- do.call(bdiag.spam, Q)
    }
    class(Q) <- "spam"
} else {
    MRA2 <- eval_basis(locs, grid, use_spam = FALSE)
    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid, use_spam = FALSE)
    MRA <- eval_basis(cbind(MRA1$W %*% alpha_x1, MRA1$W %*% alpha_y1), grid, use_spam = FALSE)

    Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=rep(0.9, length(MRA2$n_dims)), use_spam = FALSE)
    if (length(MRA2$n_dims) > 1) {
        for (m in 1:M){
            Q2[[m]] <- 2^(2*(m-1)) * Q2[[m]]
        }
        Q2 <- do.call(bdiag, Q2)
    }

    Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=rep(0.9, length(MRA1$n_dims)), use_spam = FALSE)
    if (length(MRA1$n_dims) > 1) {
        for (m in 1:M){
            Q1[[m]] <- 2^(2*(m-1)) * Q1[[m]]
        }
        Q1 <- do.call(bdiag, Q1)
    }

    Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=rep(0.9, length(MRA1$n_dims)), use_spam=FALSE)
    if (length(MRA$n_dims) > 1) {
        for (m in 1:M){
            Q[[m]] <- 2^(2*(m-1)) * Q[[m]]
        }
        Q <- do.call(bdiag, Q)
    }
}

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

i=1

##
## end setup
##

# profvis::profvis({
system.time({
    for (i in 1:4) {
        message("iteration ", i)
        N <- length(y)
        # first layer w.r.t alpha
        delta <- drop(- 1/N * (y - MRA$W %*% alpha))
        # first layer w.r.t W
        delta_W <- (2 / N * (MRA$W %*% alpha - y)) %*% t(alpha)

        # second layer
        delta_x1 <- rowSums(delta_W * (MRA$dW * MRA$ddistx))
        delta_y1 <- rowSums(delta_W * (MRA$dW * MRA$ddisty))

        # third layer
        delta2 <- delta_x1 %*% t(alpha_x1) + delta_y1 %*% t(alpha_y1) # chain rule gradient tape
        delta_x2 <- rowSums(delta2 * MRA1$dW * MRA1$ddistx)
        delta_y2 <- rowSums(delta2 * MRA1$dW * MRA1$ddisty)

        # update the gradient with adam with a penalty
        if (penalized) {
            grad <- list(t(MRA$W) %*% delta + Q %*% alpha / N,
                         t(MRA1$W) %*% delta_x1 + Q1 %*% alpha_x1 / N,
                         t(MRA1$W) %*% delta_y1 + Q1 %*% alpha_y1 / N,
                         t(MRA2$W) %*% delta_x2 + Q2 %*% alpha_x2 / N,
                         t(MRA2$W) %*% delta_y2 + Q2 %*% alpha_y2 / N)
        } else {
            # update the gradient without a penalty term
            grad <- list(t(MRA$W) %*% delta,
                         t(MRA1$W) %*% delta_x1,
                         t(MRA1$W) %*% delta_y1,
                         t(MRA2$W) %*% delta_x2,
                         t(MRA2$W) %*% delta_y2)
        }

        adam_out <- adam(i, grad, m, v)
        alpha    <- alpha    - learn_rate * adam_out$m_hat[[1]] / (sqrt(adam_out$v_hat[[1]]) + adam_out$epsilon)
        alpha_x1 <- alpha_x1 - learn_rate * adam_out$m_hat[[2]] / (sqrt(adam_out$v_hat[[2]]) + adam_out$epsilon)
        alpha_y1 <- alpha_y1 - learn_rate * adam_out$m_hat[[3]] / (sqrt(adam_out$v_hat[[3]]) + adam_out$epsilon)
        alpha_x2 <- alpha_x2 - learn_rate * adam_out$m_hat[[4]] / (sqrt(adam_out$v_hat[[4]]) + adam_out$epsilon)
        alpha_y2 <- alpha_y2 - learn_rate * adam_out$m_hat[[5]] / (sqrt(adam_out$v_hat[[5]]) + adam_out$epsilon)


        # the forward pass on the sgMRA
        MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid, use_spam = use_spam)
        MRA <- eval_basis(cbind(MRA1$W %*% alpha_x1, MRA1$W %*% alpha_y1), grid, use_spam = use_spam)
    }
})



bm  <- microbenchmark::microbenchmark(
    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid, use_spam=TRUE),
    MRA1 <- eval_basis(cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2), grid, use_spam=FALSE),
    times = 5)

bm
autoplot(bm)
##
## explore eval_basis code speed ----
##



profvis::profvis({
    for (i in 1:5){
        message("iteration ", i)
        locs = cbind(MRA2$W %*% alpha_x2, MRA2$W %*% alpha_y2)
        n_padding     = 5L
        n_neighbors   = 68
        max_points    = NULL
        basis_type    = "wendland"
        use_spam      = TRUE


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
system.time({
        for (m in 1:M) {


            # rewriten to include dW and ddist functions for more efficiency,
            # TODO: rewrite to include basis and dbasis


            D <- distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                                              radius = grid$radius[m],
                                              n_neighbors = n_neighbors)
            D$basis <- make_basis(D$V, grid$radius[m], basis_type='wendland')
            D$dbasis <- dwendland_basis(D$V, grid$radius[m])

            N_grid <- nrow(grid$locs_grid[[m]])

            if (use_spam) {
                ## use the spam sparse matrix package
                W[[m]] <- spam(D[c('ind', 'basis')], nrow =N, ncol = N_grid)
                dW[[m]] <- spam(D[c("ind", "dbasis")], nrow = N, ncol = N_grid)
                ddistx[[m]] <- spam(D[c("ind", "ddistx")], nrow = N, ncol = N_grid)
                ddisty[[m]] <- spam(D[c("ind", "ddisty")], nrow = N, ncol = N_grid)

                # }
            } else {
                ## use the Matrix sparse matrix package
                W[[m]] <- sparseMatrix(i = D$ind[, 1], j=D$ind[, 2], x = as.numeric(D$basis), dims = c(N, N_grid))
                dW[[m]] <- sparseMatrix(i = D$ind[, 1], j=D$ind[, 2], x = as.numeric(D$dbasis), dims = c(N, N_grid))
                ddistx[[m]] <- sparseMatrix(i = D$ind[, 1], j=D$ind[, 2], x = as.numeric(D$ddistx), dims = c(N, N_grid))
                ddisty[[m]] <- sparseMatrix(i = D$ind[, 1], j=D$ind[, 2], x = as.numeric(D$ddisty), dims = c(N, N_grid))

            }
        }
})

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
    }

})
