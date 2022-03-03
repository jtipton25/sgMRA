# TODO
# - work on fitting these using stochastic gradient descent or elliptical slice sampling
# - fast sparse distances: https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r
# - in python: https://stackoverflow.com/questions/56889635/is-there-a-better-and-faster-way-to-convert-from-scipy-condensed-distance-matrix
#       - python/Julia pdist: https://stackoverflow.com/questions/66237199/alternate-approach-for-pdist-from-scipy-in-julia
#       - Julia pdist
# - speedup MRA building function as this seems to be limiting the sgMRA


# Deep BayesMRA
library(spam)
library(Matrix)
library(igraph)
library(tidyverse)
library(BayesMRA)
library(patchwork)

library(sgMRA)

set.seed(42)

N <- 2^10

M <- 2
n_coarse_grid <- 10
# N <- 2^12
# M <- 1
# n_coarse_grid <- 80
# source("~/sgMRA/R/sim-deep-mra.R")
dat_sim <- sim_deep_mra(N, M, n_coarse_grid, n_layers = 3, sigma = 0.1)

y=dat_sim$y
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

p1 / p2


# make this into a function
dat <- data.frame(x = locs$x, y = locs$y,
                  layer = rep(c(1, 1, 2, 2, 3), each=N),
                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(dat_sim$MRA[[3]]$W %*% dat_sim$alpha_x[[2]], dat_sim$MRA[[3]]$W %*% dat_sim$alpha_y[[2]],
                        dat_sim$MRA[[2]]$W %*% dat_sim$alpha_x[[1]], dat_sim$MRA[[2]]$W %*% dat_sim$alpha_y[[1]],
                        dat_sim$MRA[[1]]$W %*% dat_sim$alpha))
                  # layer = rep(c(1, 1, 2), each=N),
                  # group = rep(c("x", "y", "z"), each = N),
                  # z = c(dat_sim$MRA[[2]]$W %*% dat_sim$alpha_x[[1]], dat_sim$MRA[[2]]$W %*% dat_sim$alpha_y[[1]],
                  #       dat_sim$MRA[[1]]$W %*% dat_sim$alpha))


p_layers_sim <- ggplot(dat, aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_grid(layer ~ group) +
    ggtitle("simulated layers")

p_layers_sim

# log-likelihood
1 / N * sum((dat_sim$y - dat_sim$MRA[[1]]$W %*% dat_sim$alpha)^2)

# Fit the model using sgd ----
n_iter = 2000

# add in Adam optimization schedule
out <- fit_sgd(y=dat_sim$y,
               locs=dat_sim$locs,
               grid=dat_sim$grid,
               alpha=NULL,
               alpha_x1=NULL,
               alpha_y1=NULL,
               alpha_x2=NULL,
               alpha_y2=NULL,
               learn_rate = 0.1,
               n_iter = n_iter,
               n_message = 50,
               penalized = TRUE,
               plot_during_fit = TRUE,
               use_spam=FALSE)

# re-fit the model with the current state to continue the learning

out <- fit_sgd(y=dat_sim$y,
               locs=dat_sim$locs,
               grid=dat_sim$grid,
               alpha=out$alpha,
               alpha_x1=out$alpha_x1,
               alpha_y1=out$alpha_y1,
               alpha_x2=out$alpha_x2,
               alpha_y2=out$alpha_y2,
               learn_rate = 0.01,
               n_iter = n_iter,
               n_message = 50,
               penalized = TRUE,
               plot_during_fit = TRUE,
               use_spam=FALSE,
               adam_pars = out$adam_pars)
plot(out$loss, type = 'l')

# examine the fitted layers
dat <- data.frame(x = locs$x, y = locs$y,
                  layer = rep(c(1, 1, 2, 2, 3), each=N),
                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(drop(out$MRA2$W %*% out$alpha_x2), drop(out$MRA2$W %*% out$alpha_y2),
                        drop(out$MRA1$W %*% out$alpha_x1), drop(out$MRA1$W %*% out$alpha_y1),
                        drop(out$MRA$W %*% out$alpha)))
p_layers_fit <- ggplot(dat, aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_grid(layer ~ group) +
    ggtitle("fitted layers")
p_layers_sim / p_layers_fit



dat <- data.frame(x = locs$x, y = locs$y,
                  z = z,
                  z_fit = drop(out$MRA$W %*% out$alpha))

p_fit <- ggplot(dat, aes(x = x, y = y, fill = z_fit)) +
    geom_raster() +
    scale_fill_viridis_c()
p1 / p_fit


# dat <- data.frame(z = dat_sim$z,
#                   z_fit = out$MRA$W %*% out$alpha)
p5 <- ggplot(dat, aes(x = z, y = z_fit)) +
    geom_point() +
    stat_smooth(method = 'lm')
p5


cor(dat$z, dat$z_fit)


p_resid <- ggplot(dat, aes(x = x, y = y, fill = z-z_fit)) +
    geom_raster() +
    scale_fill_viridis_c()
p_resid

dat %>%
    mutate(resid = z - z_fit) %>%
    summarize(rmse=sd(resid))


# view object memory
lapply(out,  function(x) {format(object.size(x), units="MB")})
