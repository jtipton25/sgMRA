# Laplace Approximations for log-normal likelihood

# example reference https://www.researchgate.net/publication/321347289_Anisotropic_Radial_Basis_Function_Methods_for_Continental_Size_Ice_Sheet_Simulations

set.seed(44)
library(tidyverse)
library(patchwork)
library(spam)

N <- 50^2
radius <- 0.25

locs <- as.matrix(expand_grid(x = seq(0, 1, length.out=sqrt(N)), y = seq(0, 1, length.out=sqrt(N))))
D <- fields::rdist(locs)
W <- BayesMRA::wendland_basis(D, radius)
Q <- BayesMRA::make_Q_alpha_2d(sqrt(N), phi = 0.75)
class(Q) <- "spam"
alpha <- drop(rmvnorm.prec(1, rep(0, N), Q * 50))
lambda <- W %*% alpha
dat <- data.frame(x = locs[, 1], y = locs[, 2], lambda = lambda, z = exp(lambda))
p1 <- ggplot(dat, aes(x, y, fill = lambda)) +
    geom_raster() +
    scale_fill_viridis_c()
p2 <- ggplot(dat, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c()
p1 + p2



