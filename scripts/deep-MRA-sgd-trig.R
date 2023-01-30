# TODO
# - work on fitting these using stochastic gradient descent or elliptical slice sampling


# Deep BayesMRA
library(spam)
library(Matrix)
library(igraph)
library(tidyverse)
library(BayesMRA)
library(patchwork)
library(sgMRA)

# source("~/sgMRA/R/eval_basis.R")
# Rcpp::sourceCpp("~/sgMRA/src/dist_near_cpp.cpp")
# source("~/sgMRA/R/dwendland_basis.R")

set.seed(44)
N <- 2^12
locs <- expand_grid(x=seq(0, 1, length.out=sqrt(N)),
                    y=seq(0, 1, length.out=sqrt(N)))

idx1 <- ((locs$x > 1/4) & (locs$x <= 3/4)) & ((locs$y > 1/4) & (locs$y <= 3/4))
idx2 <- ((locs$x > 1/2) & (locs$x <= 3/4)) & ((locs$y > 1/2) & (locs$y <= 3/4))
idx3 <- ((locs$x > 1/4) & (locs$x <= 1/2)) & ((locs$y > 1/4) & (locs$y <= 1/2))
z <- cos(2*pi*locs$x) * cos(2*pi*locs$y)
z[idx1] <- z[idx1] + sin(4*pi*locs$x[idx1]) * sin(4*pi*locs$y[idx1])
z[idx2] <- z[idx2] + sin(8*pi*locs$x[idx2]) * sin(8*pi*locs$y[idx2])
z[idx3] <- z[idx3] + sin(16*pi*locs$x[idx3]) * sin(16*pi*locs$y[idx3])
z <- 2*z


M <- 1
n_coarse_grid <- 200

grid <- make_grid(locs, M = M, n_coarse_grid = n_coarse_grid)
MRA <- eval_basis(locs, grid, use_spam = FALSE)
dim(MRA$W)


sigma <- 0.05
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


# sample the data

sample_percent <- 0.8
s <- sample(1:N, N * sample_percent)
y_s <- y[s]
y_oos <- y[-s]

locs_s <- locs[s, ]
locs_oos <- locs[-s, ]

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
# source("~/sgMRA/scripts/fit-deep-MRA-sgd.R")
# source("~/sgMRA/R/adam.R")
n_iter = 2000
n_steps <- 50

# linear schedule
rate_schedule <- rep(seq(500, 1, length = n_steps), each = ceiling(n_iter / n_steps))
# plot(rate_schedule)

# exponential decay schedule
# rate_schedule <- exp(rep(seq(log(0.1), log(0.001), length = n_steps), each = ceiling(n_iter / n_steps)))
plot(rate_schedule)


message("Simulated loss:", 1 / (2 * N) * sum((y - z)^2))
# profvis::profvis({
# system.time({
out <- fit_sgd(y = y_s, locs = locs_s,
               y_obs = y_obs[s],
               z = z[s],
               grid = grid,
               learn_rate = 0.001,
               rate_schedule = rate_schedule,
               n_iter = n_iter,
               n_message = 50,
               penalized = FALSE,
               plot_during_fit = TRUE,
               noisy = FALSE)

plot(out$loss, type = 'l')

n_iter = 2000
# linear schedule
rate_schedule <- rep(seq(50, 1, length = n_steps), each = ceiling(n_iter / n_steps))


# })
# resume the GD with last model fit
message("Simulated loss:", 1 / (2 * N) * sum((y - z)^2))
out <- fit_sgd(y = y_s, locs = locs_s,
               y_obs = y_obs[s],
               z = z[s],
               grid = grid,
               alpha = out$alpha,
               alpha_x1 = out$alpha_x1,
               alpha_y1 = out$alpha_y1,
               alpha_x2 = out$alpha_x2,
               alpha_y2 = out$alpha_y2,
               learn_rate = 0.001,
               rate_schedule = rate_schedule,
               n_iter = n_iter,
               n_message = 50,
               penalized = TRUE,
               plot_during_fit = TRUE,
               adam_pars = out$adam_pars,
               noisy = FALSE)

plot(out$loss, type = 'l')

# examine the fitted layers

preds <- predict_deep_MRA(out, grid, locs)
dat <- data.frame(x = locs$x, y = locs$y,
                  # layer = rep(c(1, 1, 2), each=N),
                  # group = rep(c("x", "y", "z"), each = N),
                  layer = rep(c(1, 1, 2, 2, 3), each=N),
                  group = rep(c("x", "y", "x", "y", "z"), each = N),
                  z = c(preds$z_alpha_x2, preds$z_alpha_y2,
                        preds$z_alpha_x1, preds$z_alpha_y1,
                        preds$z))

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

p1 / p_fit
# plot(out$sigma, type = 'l')
# abline(h = sigma, col = 'red')


# generate posterior mean surface
# Walpha <- matrix(0, nrow(out$W[[1]]), ncol = length(out$sigma))
# for (k in 1:length(out$sigma)) {
#     Walpha[, k] <- drop(out$W[[k]] %*% out$alpha[k, ])
# }

dat <- data.frame(x = locs$x, y = locs$y,
                  z = z,
                  z_fit = preds$z)

# p3 <- ggplot(dat, aes(x = x, y = y, fill = z_fit)) +
#     geom_raster() +
#     scale_fill_viridis_c()
# p4 <- ggplot(dat, aes(x = x, y = y, fill = z_fit - z)) +
#     geom_raster() +
#     scale_fill_viridis_c()
# p1 + p3 + p4


dat <- data.frame(z = z,
                  z_fit = preds$z)
p5 <- ggplot(dat, aes(x = z, y = z_fit)) +
    geom_point() +
    stat_smooth(method = 'lm')
p5


cor(dat$z, dat$z_fit)

# view object memory
# lapply(out,  function(x) {format(object.size(x), units="MB")})
