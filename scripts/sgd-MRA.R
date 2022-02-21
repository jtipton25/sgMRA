# TODO
# - add in a decaying learning rate
# - add in penalty term Q_alpha_tau2 with fixed penalty tau2
#     - add in adaptive penalty term tau2?
# - add in sum-to-zero constraint
#
# Do this with BayesMRA
library(tidyverse)
library(patchwork)
library(sgMRA)
library(BayesMRA)
library(spam)
library(Matrix)
library(patchwork)
set.seed(11)

N <- 400^2

## setup the spatial process
locs <- as.matrix(
    expand.grid(
        seq(0, 1, length.out = sqrt(N)),
        seq(0, 1, length.out = sqrt(N))
    )
)
# D <- fields::rdist(locs)

## fixed effects include intercept, elevation, and latitude
X <- cbind(1, rnorm(N), locs[, 2])
# X <- cbind(1, as.vector(mvnfast::rmvn(1, rep(0, N), 3 * exp(-D / 20))), locs[, 2])
p <- ncol(X)

beta <- rnorm(ncol(X))

## MRA spatio-temporal random effect
M <- 3
n_coarse <- 30

MRA    <- mra_wendland_2d(locs, M = M, n_coarse = n_coarse, use_spam = TRUE)

# MRA    <- mra_wendland_2d(locs, M = M, n_max_fine_grid = 2^8, use_spam = TRUE)
W <- as.dgCMatrix.spam(MRA$W) # using the Matrix package is faster here because we don't need the Choleksy
# W <- MRA$W
n_dims <- MRA$n_dims
dims_idx <- MRA$dims_idx

Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(0.999, length(n_dims)), prec_model = "CAR")

tau2         <- 10 * 2^(2 * (1:M - 1))
Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2)

## initialize the random effect
## set up a linear constraint so that each resolution sums to one
A_constraint <- sapply(1:M, function(i){
    tmp = rep(0, sum(n_dims))
    tmp[dims_idx == i] <- 1
    return(tmp)
})

a_constraint <- rep(0, M)
alpha   <- as.vector(rmvnorm.prec.const(n = 1, mu = rep(0, sum(n_dims)), Q = Q_alpha_tau2, A = t(A_constraint), a = a_constraint))

sigma2 <- runif(1, 0.25, 0.5)

y <- as.numeric(X %*% beta + W %*% alpha + rnorm(N, 0, sqrt(sigma2)))

dat <- data.frame(x=locs[, 1], y = locs[, 2], z = y)

p_sim <- ggplot(dat, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Simulated")
p_sim

# gradient descent function for MRA using minibatch ----

U <- cbind(X, W)


profvis::profvis(
# system.time(
    dat_mini <- regression_gradient_descent(c(y),
                                       X, W, inits = rnorm(ncol(U)),
                                       threshold = 0.000001,
                                       alpha = 0.5, num_iters = 500,
                                       print_every = 10,
                                       minibatch_size = 2^6)
)

# Using full gradient
# profvis::profvis(
system.time(
    dat_full <- regression_gradient_descent(c(y),
                                       X, W, inits = rnorm(ncol(U)),
                                       threshold = 0.000001,
                                       alpha = 0.5, num_iters = 500,
                                       print_every = 10,
                                       minibatch_size = NULL)
)


# predictions from the minibatch model
beta_idx <- grepl("beta", names(dat_mini))
beta_fit_mini <- dat_mini[nrow(dat_mini), beta_idx]
y_pred_mini <- U %*% as.numeric(beta_fit_mini)
# predictions from the full model
beta_idx <- grepl("beta", names(dat_full))
beta_fit_full <- dat_full[nrow(dat_full), beta_idx]
y_pred_full <- U %*% as.numeric(beta_fit_full)


dat <- data.frame(x = locs[, 1], y = locs[, 2], z = drop(W %*% alpha), y_sim = y,
                  y_pred_mini = drop(y_pred_mini), y_pred_full = drop(y_pred_full))

p_mini <- ggplot(dat, aes(x, y, fill = y_pred_mini)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Minibatch fit")
p_full <- ggplot(dat, aes(x, y, fill = y_pred_full)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Full fit")

p_sim + p_mini + p_full

sd(y - y_pred_mini)
sd(y - y_pred_full)


dat %>%
    pivot_longer(cols= y_pred_mini:y_pred_full) %>%
    ggplot(aes(x = y_sim, y=value)) +
    geom_point() +
    facet_wrap(~ name)

tU <- t(U)
tUU <- t(U) %*% U
theta_test <- rnorm(ncol(U))

bm <- microbenchmark::microbenchmark(
    tUU %*% theta_test,
    tU %*% (U %*% theta_test), times = 10
)
autoplot(bm)


tUy <- t(U) %*% y
idx <- sample(1:length(y), 2^6)

y_idx <- y[idx]
tUy_idx <- tU[, idx] %*% y_idx


bm <- microbenchmark::microbenchmark(
    gradient_fun(y, tUy, tUU, theta_test),
    gradient_fun(y[idx], tU[, idx] %*% y[idx], tUU, theta_test), times = 10)
autoplot(bm)
# all.equal(    gradient_fun(y, U, theta_test),
#               gradient_fun2(y, U, tUy, tUU, theta_test))


ggplot(dat, aes(x = iteration, y = loss)) +
    geom_point() +
    geom_line() +
    ggtitle("loss as a function of iteration")


# Plot the regression parameters

# convert beta into a data.frame
# results from lm()

dat %>%
    pivot_longer(cols = 1:10, names_to = "parameter", values_to = "beta") %>%
    ggplot(aes(x = iteration, y = beta, color = parameter)) +
    geom_line() +
    ggtitle("Parameter estimates by iteration")


# layout(matrix(1:2, 2, 1))
# plot(y, cbind(X, W) %*% as.matrix(dat)[nrow(dat), 1:(ncol(dat) - 2)])
# plot(y, X %*% beta + W %*% alpha)

# fitted MSE and R^2
preds <- as.numeric(cbind(X, W) %*% unlist(dat[nrow(dat), ][1:(ncol(dat) - 2)]))
mean((y - preds)^2)
cor(y, preds)^2
# oracle MSE and R^2
mu <- as.numeric(X %*% beta + W %*% alpha)
mean((y - mu)^2)
cor(y, mu)^2

# plot spatial random field
dat <- data.frame(x = locs[, 1], y = locs[, 2], value = y, mu = mu, pred = preds)
plot_lims <- range(c(y, mu, preds))
p_sim <- ggplot(dat) +
    geom_raster(aes(x, y, fill = value)) +
    scale_fill_viridis_c(end = 0.8, limits = plot_lims) +
    ggtitle("simulated data")
p_latent <- ggplot(dat) +
    geom_raster(aes(x, y, fill = mu)) +
    scale_fill_viridis_c(end = 0.8, limits = plot_lims) +
    ggtitle("simulated latent process")
p_fit <- ggplot(dat) +
    geom_raster(aes(x, y, fill = preds)) +
    scale_fill_viridis_c(end = 0.8, limits = plot_lims) +
    ggtitle("predicted latent process")
p_sim + p_latent + p_fit
