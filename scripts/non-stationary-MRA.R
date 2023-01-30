# example reference https://www.researchgate.net/publication/321347289_Anisotropic_Radial_Basis_Function_Methods_for_Continental_Size_Ice_Sheet_Simulations

set.seed(44)
library(tidyverse)
library(sgMRA)
N_grid <- 40^2
a <- 0.5
b <- 2
radius <- 0.1
trunc <- 0.1

# plot a kernel
grid <- as.matrix(expand_grid(x = seq(0, 1, length.out=sqrt(N_grid)), y = seq(0, 1, length.out=sqrt(N_grid))))
point <- data.frame(x=0.15, y=0.5)
D <- rep(0, N_grid)
for (i in 1:N_grid) {
    # D[i] <- sqrt(sum(c(a, b) * (point - grid[i, ])^2))
    D[i] <- sqrt(sum((point - grid[i, ])^2))
}
W <- BayesMRA::wendland_basis(D, radius)
dat <- data.frame(x=grid[, 1], y=grid[, 2], W = drop(W))
ggplot(dat, aes(x=x, y=y, fill=W)) +
    geom_raster() +
    scale_fill_viridis_c()


# construct the deviations (this will be a 3-d array in the real code) ----
N <- 100^2
# N=1
# locs <- matrix(c(0.15, 0.5), 1, 2)
locs <- as.matrix(expand_grid(x = seq(0, 1, length.out=sqrt(N)), y = seq(0, 1, length.out=sqrt(N))))
sq_devs <- array(0, dim=c(N, N_grid, 2))
for (i in 1:N) {
    # can make this parallelizable
    # sq_devs[i, , ] <- (locs[i, ] - grid)^2
    sq_devs[i, , ] <- sweep(grid, 2, FUN='-', locs[i, ])^2
}

# create constant anisotropic kernels over space ----
# D <- matrix(0, N, N_grid)
# for (i in 1:N) {
#     # D[i, ] <- sqrt(rowSums(sq_devs[i, , ]))# * cbind(rep(a N), rep(b, N))))
#     D[i, ] <- sqrt(rowSums(sq_devs[i, , ] * cbind(rep(a, N_grid), rep(b, N_grid))))
# }

DD <- distance_near_nonstationary(sq_devs,
                                  rep(a, N),
                                  rep(b, N),
                                  radius = radius)
# truncate to a sparse matrix
D <- sparseMatrix(i = DD$ind[, 1], j = DD$ind[, 2], x = as.vector(DD$V), dims = dim(sq_devs)[1:2], repr = "T")


W <- D
W@x <- BayesMRA::wendland_basis(D@x, radius)
# W <- BayesMRA::wendland_basis(D, radius)


# pick four points
idx <- c(512, 1024, 1536, 2048)
dat <- data.frame(x=grid[, 1], y=grid[, 2], W = c(t(as.matrix(W)[idx, ])), location = rep(idx, each = N_grid))
ggplot(dat, aes(x=x, y=y, fill=W)) +
    geom_raster() +
    facet_wrap(~ location) +
    scale_fill_viridis_c()



# create spatially varying anisotropic kernels over space ----
# give this a much smaller grid than the overall grid
library(spam)
N_grid_kernels <- 40^2
grid_kernels <- as.matrix(expand_grid(x = seq(0, 1, length.out=sqrt(N_grid_kernels)),
                                      y = seq(0, 1, length.out=sqrt(N_grid_kernels))))

D_kernels <- matrix(0, N, N_grid_kernels)
for (i in 1:N) {
    D_kernels[i, ] <- sqrt(rowSums(sweep(grid_kernels, 2, FUN='-', locs[i, ])^2))
}

W_kernels <- BayesMRA::wendland_basis(D_kernels, radius)
Q_kernels <- sgMRA::make_Q(sqrt(N_grid_kernels), 0.99)
alpha_a <- drop(spam::rmvnorm.prec(1, rep(0, N_grid_kernels), Q_kernels)) * 0.5
alpha_b <- drop(spam::rmvnorm.prec(1, rep(0, N_grid_kernels), Q_kernels)) * 0.5
a <- pmax(drop(W_kernels %*% alpha_a + 0.5), trunc)
b <- pmax(drop(W_kernels %*% alpha_b + 0.5), trunc)


dat_sim <- data.frame(x = locs[, 1],
                      y = locs[, 2],
                      z = as.vector(W %*% alpha),
                      a = a,
                      b = b)

p_sim_a <- ggplot(dat_sim, aes(x, y, fill = a)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("simulated a")
p_sim_b <- ggplot(dat_sim, aes(x, y, fill = b)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("simulated b")

p_sim_a + p_sim_b


# D <- matrix(0, N, N_grid)
# for (i in 1:N) {
#     D[i, ] <- sqrt(rowSums(sq_devs[i, , ] * cbind(rep(a[i], N_grid), rep(b[i], N_grid))))
# }
#
#
# W <- BayesMRA::wendland_basis(D, radius)

DD <- distance_near_nonstationary(sq_devs,
                                  a,
                                  b,
                                  radius = radius)
# truncate to a sparse matrix
D <- sparseMatrix(i = DD$ind[, 1], j = DD$ind[, 2], x = as.vector(DD$V), dims = dim(sq_devs)[1:2], repr = "T")


W <- D
W@x <- BayesMRA::wendland_basis(D@x, radius)


# pick four points
idx <- c(10, 200, 512, 600, 800, 1024, 1200, 1400, 1536, 1650, 1850, 2048)
idx <- sample(1:N, 25)
dat <- data.frame(x=grid[, 1], y=grid[, 2], W = c(t(as.matrix(W)[idx, ])), location = rep(idx, each = N_grid))
ggplot(dat, aes(x=x, y=y, fill=W)) +
    geom_raster() +
    facet_wrap(~ location) +
    scale_fill_viridis_c()

# Plot a realization of the simulated process ----
Q <- sgMRA::make_Q(sqrt(N_grid), 0.99)
alpha <- spam::rmvnorm.prec(16, rep(0, N_grid), Q) * 0.1
dat <- data.frame(x=locs[, 1], y=locs[, 2], Z = as.vector(W %*% t(alpha)), replicate = rep(1:16, each=N))
ggplot(dat, aes(x=x, y=y, fill=Z)) +
    geom_raster() +
    facet_wrap(~ replicate) +
    scale_fill_viridis_c()


alpha <- drop(alpha[1, ])
z <- as.vector(W %*% alpha)
epsilon <- rnorm(N, 0, 0.1)
y <- z + epsilon


# Fit the nonstationary process

n_iter <- 500
n_steps <- 50
n_message = 10

rate_schedule <- rep(seq(0.05, 0.001, length = n_steps), each = ceiling(n_iter / n_steps))
# plot(rate_schedule)

# exponential decay schedule
# rate_schedule <- exp(rep(seq(log(0.1), log(0.001), length = n_steps), each = ceiling(n_iter / n_steps)))
plot(rate_schedule)



message("oracle loss = ", 1 / (2 * N) * sum((y - W %*% alpha)^2))
# add in noisy gradients
# need to add in sparse matrices
# figure out how to update the radius and trunc parameters as well
out <- fit_nonstationary_MRA(y, locs, grid, grid_kernels, n_iter=n_iter, rate_schedule = rate_schedule,
                      n_message = n_message, plot_during_fit = TRUE, radius = radius, trunc = trunc)


# profvis::profvis(fit_nonstationary_MRA(y, locs, grid, grid_kernels, n_iter=n_iter, learn_rate = rate_schedule,
#                       n_message = n_message, plot_during_fit = TRUE))

n_iter <- 1000
n_steps <- 20
n_message = 10

rate_schedule <- rep(seq(0.5, 0.001, length = n_steps), each = ceiling(n_iter / n_steps))
out <- fit_nonstationary_MRA(y, locs, grid, grid_kernels, n_iter=n_iter, rate_schedule = rate_schedule,
                             n_message = n_message, plot_during_fit = TRUE,
                             alpha = out$alpha,
                             alpha_a = out$alpha_a,
                             alpha_b = out$alpha_b,
                             adam_pars = out$adam_pars,
                             radius = radius, trunc = trunc)



# plot the fitted a and b vs. the simulated a and b
dat_sim <- data.frame(x = locs[, 1],
                      y = locs[, 2],
                      z = as.vector(W %*% alpha),
                      a = a,
                      b = b)
dat_fitted <- data.frame(x = locs[, 1],
                         y = locs[, 2],
                         z = as.vector(out$W %*% out$alpha),
                         a = as.vector(pmax(out$W_kernels %*% out$alpha_a, trunc)),
                         b = as.vector(pmax(out$W_kernels %*% out$alpha_b, trunc)))


p_sim_a <- ggplot(dat_sim, aes(x, y, fill = a)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("simulated a")
p_fitted_a <- ggplot(dat_fitted, aes(x, y, fill = a)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("fitted a")

p_sim_b <- ggplot(dat_sim, aes(x, y, fill = b)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("simulated b")
p_fitted_b <- ggplot(dat_fitted, aes(x, y, fill = b)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("fitted b")

(p_sim_a + p_fitted_a) / (p_sim_b + p_fitted_b)


p_sim_z <- ggplot(dat_sim, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("simulated z")
p_fitted_z <- ggplot(dat_fitted, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("fitted z")

p_sim_z / p_fitted_z


dat_all <-
    dat_sim |>
    mutate(type = "sim") |>
    rbind(mutate(dat_fitted, type = 'fitted')) |>
    pivot_wider(names_from = type, names_glue = "{type}_{.value}", values_from = c("z", "a", "b")) |>
    as_tibble()

ggplot(dat_all) +
    geom_point(aes(x = fitted_z, y = sim_z))

ggplot(dat_all) +
    geom_raster(aes(x, y, fill = fitted_z - sim_z)) +
    scale_fill_viridis_c() +
    ggtitle("resitual z")

# Fit using MCMC ----
# source(here::here("R", "nonstationary_MRA_ess.R"))

# out <- nonstationary_MRA_ess(y, locs, grid, grid_kernels, radius=radius, trunc = 0.01, n_mcmc = 10, n_message = 5)
