# TODO
# - work on fitting these using stochastic gradient descent or elliptical slice sampling


# Deep BayesMRA
library(spam)
library(Matrix)
library(igraph)
library(tidyverse)
library(BayesMRA)
library(patchwork)

source("~/sgMRA/R/eval_basis.R")
source("~/sgMRA/R/dwendland_basis.R")

set.seed(44)
N <- 100^2
locs <- expand_grid(x=seq(0, 1, length.out=sqrt(N)),
                    y=seq(0, 1, length.out=sqrt(N)))



M <- 1
n_coarse_grid <- 40

grid <- make_grid(locs, M = M, n_coarse_grid = n_coarse_grid)

# construct the first layer
# MRA1 <- mra_wendland_2d(as.matrix(locs), M = M, n_coarse_grid = n_coarse_grid)
MRA1 <- eval_basis(locs, grid)
W1 <- MRA1$W
Q1 <- make_Q_alpha_2d(sqrt(MRA1$n_dims), phi=0.9)
class(Q1) <- "spam"
alpha_x1 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=rep(1, nrow(W1)) %*% W1, a=0))
alpha_y1 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=rep(1, nrow(W1)) %*% W1, a=0))
# alpha_x1 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=t(rep(1, ncol(W1))), a=0))
# alpha_y1 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q1)), Q = Q1, A=t(rep(1, ncol(W1))), a=0))

# construct the second layer
# MRA2 <- mra_wendland_2d(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), M = M, n_coarse_grid = n_coarse_grid)
MRA2 <- eval_basis(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), grid)
# MRA2 <- mra_wendland_2d_pred(cbind(W1 %*% alpha_x1, W1 %*% alpha_y1), MRA1)
W2 <- MRA2$W
Q2 <- make_Q_alpha_2d(sqrt(MRA2$n_dims), phi=0.9)
class(Q2) <- "spam"
alpha_x2 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=rep(1, nrow(W2)) %*% W2, a=0))
alpha_y2 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=rep(1, nrow(W2)) %*% W2, a=0))
# alpha_x2 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=t(rep(1, ncol(W2))), a=0))
# alpha_y2 <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q2)), Q = Q2, A=t(rep(1, ncol(W2))), a=0))

# construct the final layer
# MRA <- mra_wendland_2d(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), M = M, n_coarse_grid = n_coarse_grid)
MRA <- eval_basis(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), grid)
# MRA <- mra_wendland_2d_pred(cbind(W2 %*% alpha_x2, W2 %*% alpha_y2), MRA1)
W <- MRA$W
Q <- make_Q_alpha_2d(sqrt(MRA$n_dims), phi=0.9)
class(Q) <- "spam"
alpha <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q)), Q = Q, A=rep(1, nrow(W)) %*% W, a=0))
# alpha <- drop(rmvnorm.prec.const(1, mu = rep(0, nrow(Q)), Q = Q, A=t(rep(1, ncol(W))), a=0))
z <- W %*% alpha

sigma <- 0.5
epsilon <- rnorm(N, 0, sigma)
y_obs <- z + epsilon
y <- y_obs
dat <- data.frame(x = locs$x, y = locs$y, z = z, y_obs = y_obs)


p1 <- ggplot(dat, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c()

p2 <- ggplot(dat, aes(x = x, y = y, fill = y_obs)) +
    geom_raster() +
    scale_fill_viridis_c()

p1 + p2


dat <- data.frame(x = locs$x, y = locs$y, layer = rep(c(1, 1, 2, 2, 3), each=N),
                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(W1 %*% alpha_x1, W1 %*% alpha_y1,
                        W2 %*% alpha_x2, W2 %*% alpha_y2,
                        W %*% alpha))
ggplot(dat, aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_grid(layer~ group)


# Fit the MCMC using ESS ----
source("~/sgMRA/scripts/deep_mra_ess.R")
source("~/sgMRA/R/ess.R")
out <- deep_mra_ess(y_obs, locs,
                    alpha_x1 = alpha_x1,
                    sample_alpha_x1 = FALSE,
                    alpha_y1 = alpha_y1,
                    sample_alpha_y1 = FALSE,
                    alpha_x2 = alpha_x2,
                    sample_alpha_x2 = FALSE,
                    # alpha_y2 = alpha_y2,
                    alpha_y2 = NULL,
                    sample_alpha_y2 = TRUE,
                    alpha = alpha,
                    # alpha = NULL,
                    sample_alpha = FALSE,
                    n_message = 1, n_mcmc=500,
                    M=M, n_coarse_grid = n_coarse_grid)

# source("~/sgMRA/scripts/deep_mra_mh.R")
# source("~/sgMRA/R/update-tuning.R")
# Fit the MCMC using MH
out <- deep_mra_mh(y_obs, locs,
                   alpha_x1 = alpha_x1,
                   sample_alpha_x1 = FALSE,
                   alpha_y1 = alpha_y1,
                   sample_alpha_y1 = FALSE,
                   alpha_x2 = alpha_x2,
                   sample_alpha_x2 = FALSE,
                   # alpha_y2 = alpha_y2,
                   alpha_y2 = NULL,
                   sample_alpha_y2 = TRUE,
                   alpha = alpha,
                   # alpha = NULL,
                   sample_alpha = FALSE,
                    n_message = 10, n_mcmc = 500,
                   M=M, n_coarse_grid = n_coarse_grid)


plot(out$sigma, type = 'l')
abline(h = sigma, col = 'red')


# generate posterior mean surface
# Walpha <- matrix(0, nrow(out$W[[1]]), ncol = length(out$sigma))
# for (k in 1:length(out$sigma)) {
#     Walpha[, k] <- drop(out$W[[k]] %*% out$alpha[k, ])
# }

dat <- data.frame(x = locs$x, y = locs$y,
                  z = z,
                  z_mean = apply(out$Walpha[(length(out$sigma)/2):length(out$sigma), ], 2, mean),
                  z_last = out$Walpha[length(out$sigma),])

p3 <- ggplot(dat, aes(x = x, y = y, fill = z_mean)) +
    geom_raster() +
    scale_fill_viridis_c()
p33 <- ggplot(dat, aes(x = x, y = y, fill = z_last)) +
    geom_raster() +
    scale_fill_viridis_c()
p333 <- ggplot(dat, aes(x = x, y = y, fill = z_mean - z)) +
    geom_raster() +
    scale_fill_viridis_c()
p1 + p3


dat <- data.frame(z = z,
                  z_mean = apply(out$Walpha[(length(out$sigma)/2):length(out$sigma), ], 2, mean),
                  z_last = out$Walpha[length(out$sigma),])
p4 <- ggplot(dat, aes(x = z, y = z_mean)) +
    geom_point() +
    stat_smooth(method = 'lm')
p5 <- ggplot(dat, aes(x = z, y = z_last)) +
    geom_point() +
    stat_smooth(method = 'lm')
p4 + p5


cor(dat$z, dat$z_mean)
cor(dat$z, dat$z_last)

# view object memory
lapply(out,  function(x) {format(object.size(x), units="MB")})
