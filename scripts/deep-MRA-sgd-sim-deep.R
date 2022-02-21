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
N <- 2^8
locs <- expand_grid(x=seq(0, 1, length.out=sqrt(N)),
                    y=seq(0, 1, length.out=sqrt(N)))

idx1 <- ((locs$x > 1/4) & (locs$x <= 3/4)) & ((locs$y > 1/4) & (locs$y <= 3/4))
idx2 <- ((locs$x > 1/2) & (locs$x <= 3/4)) & ((locs$y > 1/2) & (locs$y <= 3/4))
idx3 <- ((locs$x > 1/4) & (locs$x <= 1/2)) & ((locs$y > 1/4) & (locs$y <= 1/2))
z <- cos(2*pi*locs$x) * cos(2*pi*locs$y)
z[idx1] <- z[idx1] + sin(4*pi*locs$x[idx1]) * sin(4*pi*locs$y[idx1])
z[idx2] <- z[idx2] + sin(8*pi*locs$x[idx2]) * sin(8*pi*locs$y[idx2])
z[idx3] <- z[idx3] + sin(16*pi*locs$x[idx3]) * sin(16*pi*locs$y[idx3])



M <- 1
n_coarse_grid <- 20

grid <- make_grid(locs, M = M, n_coarse_grid = n_coarse_grid)

sigma <- 0.02
epsilon <- rnorm(N, 0, sigma)
y_obs <- z + epsilon
y <- y_obs
dat <- data.frame(x = locs$x, y = locs$y, z = z, y_obs = y_obs)


p1 <- ggplot(dat, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    # scale_fill_distiller(palette = "RdYlBu")
    scale_fill_viridis_c()
p2 <- ggplot(dat, aes(x = x, y = y, fill = y_obs)) +
    geom_raster() +
    # scale_fill_distiller(palette = "RdYlBu")
    scale_fill_viridis_c()

p1 + p2


# dat <- data.frame(x = locs$x, y = locs$y,
#                   # layer = rep(c(1, 1, 2, 2, 3), each=N),
#                   layer = rep(c(1, 1, 2), each=N),
#                   group = rep(c("x", "y", "z"), each = N),
#                   # group = rep(c("x", "y", "x", "y", "z"), each = N),
#                   z = c(W1 %*% alpha_x1, W1 %*% alpha_y1,
#                         # W2 %*% alpha_x2, W2 %*% alpha_y2,
#                         W %*% alpha))
# p_layers_sim <- ggplot(dat, aes(x, y, fill=z)) +
#     geom_raster() +
#     scale_fill_viridis_c() +
#     facet_grid(layer ~ group) +
#     ggtitle("simulated layers")
#
# p_layers_sim

# Fit the model using sgd ----
source("~/sgMRA/scripts/fit-deep-MRA-sgd.R")
source("~/sgMRA/R/adam.R")
n_iter = 10000
rate_schedule = rep(seq(0.25, 0.01, length = 100), each = ceiling(n_iter / 100))

# add in Adam optimization schedule
# profvis::profvis(
system.time(
    out <- fit_sgd(y = y, locs = locs, grid = grid,
               alpha=NULL,
               # alpha=alpha,
               alpha_x1=NULL,
               # alpha_x1=NULL,
               alpha_y1=NULL,
               # alpha_y1=NULL,
               # alpha_x2=NULL,
               # # alpha_x2=NULL,
               # alpha_y2=NULL,
               # # alpha_y2=NULL,
               learn_rate = 0.01,
               # rate_schedule = rate_schedule,
               n_iter = n_iter,
               n_message = 100)
)
# M=1, n_coarse_grid=10, layers = 3, fit-time:  elapsed, loss:
# M=1, n_coarse_grid=10, layers = 2, fit-time: 364.731 elapsed, loss: 0.32
# M=2, n_coarse_grid=30, layers = 2, fit-time: 1178.980 elapsed, loss: 0.0741
# )
# plot(out$alpha, alpha)
# abline(0, 1, col = 'red')
plot(out$loss, type = 'l')

# examine the fitted layers
dat <- data.frame(x = locs$x, y = locs$y,
                  # layer = rep(c(1, 1, 2, 2, 3), each=N),
                  layer = rep(c(1, 1, 2), each=N),
                  group = rep(c("x", "y", "z"), each = N),
                  # group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(out$MRA1$W %*% out$alpha_x1, out$MRA1$W %*% out$alpha_y1,
                        # out$MRA2$W %*% out$alpha_x2, out$MRA2$W %*% out$alpha_y2,
                        out$MRA$W %*% out$alpha))
# plot_func <- function(df, name) {
#     ggplot(data = df, aes(x = x, y = y, fill = z)) +
#         geom_raster() +
#         scale_fill_viridis_c() +
#         scale_fill_continuous(name = name)
# }
#
# nested_tmp <- dat %>%
#     group_by(group, layer) %>%
#     nest() %>%
#     mutate(plots = map2(data, group, plot_func)) %>%
#     gridExtra::grid.arrange(grobs = .$plots)


p_layers_fit <- ggplot(dat, aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_grid(layer ~ group) +
    ggtitle("fitted layers")
p1 / p_layers_fit

p_fit <- dat %>%
    filter(group == "z") %>%
    ggplot(aes(x, y, fill=z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("fitted process")

p1 + p_fit
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

p3 <- ggplot(dat, aes(x = x, y = y, fill = z_fit)) +
    geom_raster() +
    scale_fill_viridis_c()
p4 <- ggplot(dat, aes(x = x, y = y, fill = z_fit - z)) +
    geom_raster() +
    scale_fill_viridis_c()
p1 + p3 + p4


dat <- data.frame(z = z,
                  z_fit = out$MRA$W %*% out$alpha)
p5 <- ggplot(dat, aes(x = z, y = z_fit)) +
    geom_point() +
    stat_smooth(method = 'lm')
p5


cor(dat$z, dat$z_fit)

# view object memory
lapply(out,  function(x) {format(object.size(x), units="MB")})
