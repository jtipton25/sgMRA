# TODO
# - add in a decaying learning rate
# - add in penalty term Q_alpha_tau2 with fixed penalty tau2
#     - add in adaptive penalty term tau2?
# - add in sum-to-zero constraint
#
library(tidyverse)
library(patchwork)
library(sgMRA)
set.seed(404)


dat_sim <- sim_example_splines()

# plot the simulated data example ----
layout(matrix(1:2, 2, 1))
plot(dat_sim$X[, 2], dat_sim$y - dat_sim$Z %*% dat_sim$alpha)
abline(dat_sim$beta, col = 'red')
plot(seq(0, 1, length.out = length(dat_sim$y)), dat_sim$y)
lines(seq(0, 1, length.out = length(dat_sim$y)),
      dat_sim$Z %*% dat_sim$alpha + dat_sim$X[, 1] * dat_sim$beta[1],
      col = "red", lwd = 2)






# gradient descent function for regression ----

dat <- regression_gradient_descent(y=dat_sim$y, X=dat_sim$X, W=dat_sim$Z,
                                   threshold = 0.0000001, alpha = 0.1,
                                   num_iters = 50000, print_every = 5000,
                                   # minibatch_size = 2^6
                                   minibatch_size = NULL)

dat %>%
    filter(row_number() %% 100 == 0) %>% # if needing to reduce samples
    ggplot(aes(x = iteration, y = loss)) +
    geom_point() +
    geom_line() +
    ggtitle("loss as a function of iteration")


# Fit the model with lm() ----
model <- lm(dat_sim$y ~ cbind(dat_sim$X, dat_sim$Z)-1)
dat_coef <- data.frame(beta = coef(model), parameter = paste0("beta[", 1:(ncol(dat_sim$X) + ncol(dat_sim$Z)), "]"))


# Plot the regression parameters ----
dat %>%
    filter(row_number() %% 100 == 0) %>% # if needing to reduce samples
    pivot_longer(cols = 1:10, names_to = "parameter", values_to = "beta") %>%
    ggplot(aes(x = iteration, y = beta, color = parameter)) +
    geom_line() +
    geom_hline(data = dat_coef, aes(yintercept = beta, color = parameter), lty = 2) +
    theme(legend.position = "none") +
    ggtitle("Parameter estimates by iteration")




model$coefficients
tail(dat, n=1)[1:10]


# fitted MSE and R^2
preds <- as.numeric(cbind(dat_sim$X, dat_sim$Z) %*% unlist(dat[nrow(dat), ][1:(ncol(dat) - 2)]))
mean((dat_sim$y - preds)^2)
cor(dat_sim$y, preds)^2
# oracle MSE and R^2
mu <- as.numeric(dat_sim$X %*% dat_sim$beta + dat_sim$Z %*% dat_sim$alpha)
mean((dat_sim$y - mu)^2)
cor(dat_sim$y, mu)^2

ggplot(data.frame(y = dat_sim$y, y_pred = preds)) +
    geom_point(aes(y_pred, y)) +
    geom_abline(color = 'red')

# p <- ggplot(data = NULL, aes(x = X[, 2], y = y - Z %*% alpha)) +
#     geom_point() +
#     geom_abline(data = dat, aes(intercept = `beta[0]`, slope = `beta[1]`, color = iteration), alpha = 0.25) +
#     # geom_abline(intercept = dat$`beta[0]`, slope = dat$`beta[1]`, color = dat$iteration) +
#     scale_color_viridis_c(direction = -1)
# p

