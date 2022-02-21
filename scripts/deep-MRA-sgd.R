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

N <- 2^12
M <- 1
n_coarse_grid <- 80
source("~/sgMRA/R/sim-deep-mra.R")
dat_sim <- sim_deep_mra(N, M, n_coarse_grid, n_layers = 3, sigma = 0.01)
str(dat_sim)
length(dat_sim$MRA)
str(dat_sim$alpha)
str(dat_sim$alpha_x)

locs <- dat_sim$locs
z <- dat_sim$z
y_obs<- dat_sim$y

dat <- data.frame(x = locs$x, y = locs$y, z = z, y_obs = y_obs)


p1 <- ggplot(dat, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c()

p2 <- ggplot(dat, aes(x = x, y = y, fill = y_obs)) +
    geom_raster() +
    scale_fill_viridis_c()

p1 + p2


# make this into a function
dat <- data.frame(x = locs$x, y = locs$y,
                  layer = rep(c(1, 1, 2, 2, 3), each=N),
                  # layer = rep(c(1, 1, 2), each=N),
                  # group = rep(c("x", "y", "z"), each = N),
                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(dat_sim$MRA[[3]]$W %*% dat_sim$alpha_x[[2]], dat_sim$MRA[[3]]$W %*% dat_sim$alpha_y[[2]],
                        dat_sim$MRA[[2]]$W %*% dat_sim$alpha_x[[1]], dat_sim$MRA[[2]]$W %*% dat_sim$alpha_y[[1]],
                        dat_sim$MRA[[1]]$W %*% dat_sim$alpha))


p_layers_sim <- ggplot(dat, aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_grid(layer ~ group) +
    ggtitle("simulated layers")

p_layers_sim

# Fit the model using sgd ----
source("~/sgMRA/scripts/fit-deep-MRA-sgd.R")
source("~/sgMRA/R/adam.R")
n_iter = 1000

# add in Adam optimization schedule
out <- fit_sgd(y=dat_sim$y,
               locs=dat_sim$locs,
               grid=dat_sim$grid,
               # alpha=alpha,
               alpha=NULL,
               alpha_x1=NULL,
               # alpha_x1=NULL,
               alpha_y1=NULL,
               # alpha_y1=NULL,
               alpha_x2=NULL,
               alpha_y2=NULL,
               learn_rate = 0.001,
               n_iter = n_iter,
               n_message = 1)

plot(out$loss, type = 'l')

# plot(out$alpha, alpha)
# abline(0, 1, col = 'red')

# plot first layer
# plot(out$MRA$W %*% out$alpha, W %*% alpha)
# abline(0, 1, col = 'red')

# plot second layer x component
# plot(out$MRA1$W %*% out$alpha_x1, W1 %*% alpha_x1)
# abline(0, 1, col = 'red')

# plot second layer y component
# plot(out$MRA1$W %*% out$alpha_y1, W1 %*% alpha_y1)
# abline(0, 1, col = 'red')

# takes too long to plot

# examine the fitted layers
dat <- data.frame(x = locs$x, y = locs$y,
                  layer = rep(c(1, 1, 2, 2, 3), each=N),
                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(out$MRA1$W %*% out$alpha_x1, out$MRA1$W %*% out$alpha_y1,
                        out$MRA2$W %*% out$alpha_x2, out$MRA2$W %*% out$alpha_y2,
                        out$MRA$W %*% out$alpha))
p_layers_fit <- ggplot(dat, aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_grid(layer ~ group) +
    ggtitle("fitted layers")
p_layers_sim / p_layers_fit
# plot(out$sigma, type = 'l')
# abline(h = sigma, col = 'red')


# generate posterior mean surface
# Walpha <- matrix(0, nrow(out$W[[1]]), ncol = length(out$sigma))
# for (k in 1:length(out$sigma)) {
#     Walpha[, k] <- drop(out$W[[k]] %*% out$alpha[k, ])
# }

dat <- data.frame(x = locs$x, y = locs$y,
                  z = z,
                  z_fit = out$MRA$W %*% out$alpha)

p_fit <- ggplot(dat, aes(x = x, y = y, fill = z_fit)) +
    geom_raster() +
    scale_fill_viridis_c()
p1 + p_fit


dat <- data.frame(z = dat_sim$z,
                  z_fit = out$MRA$W %*% out$alpha)
p5 <- ggplot(dat, aes(x = z, y = z_fit)) +
    geom_point() +
    stat_smooth(method = 'lm')
p5


cor(dat$z, dat$z_fit)

# view object memory
lapply(out,  function(x) {format(object.size(x), units="MB")})