library(sgMRA)
library(purrr)
library(furrr)
# Deep BayesMRA
library(spam)
library(Matrix)
library(igraph)
library(tidyverse)
library(BayesMRA)
library(patchwork)

library(sgMRA)

set.seed(404)

N <- 2^16

M <- 3
n_coarse_grid <- 80
# N <- 2^12
# M <- 1
# n_coarse_grid <- 80
# source("~/sgMRA/R/sim-deep-mra.R")
dat_sim <- sim_deep_mra(N, M, n_coarse_grid, n_layers = 3, sigma = 0.1, use_spam = FALSE)

# fun <- function(ncores) {
#     sim_deep_mra(N, M, n_coarse_grid, n_layers = 3, sigma = 0.1, use_spam = FALSE, ncores = ncores)
# }
# bm <- microbenchmark::microbenchmark(
#     fun(1L),
#     fun(4L),
#     fun(16L),
#     fun(32L),
#     times = 4)
# bm
# print(bm, unit='relative')
# autoplot(bm)

# profvis::profvis({
#     sim_deep_mra(N, M, n_coarse_grid, n_layers = 3, sigma = 0.1, use_spam = FALSE, ncores = 32L)
# })


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


m=3
radius <- 1


calc_dist_par <- function(i, m, locs, grid) {
    tmp=as.data.frame(t(distance_near_row_cpp(i=as.integer(i-1), locs = as.matrix(locs[i,]), locs_grid = as.matrix(grid$locs_grid[[m]]), radius=grid$radius[m])))
    colnames(tmp) = c("I", "J", "V", "ddistx", "ddisty")
    tmp
}

tst=calc_dist_par(2, 1, locs, grid)


# 20 secs
D <- vector(mode = 'list', length =  M)
# for (m in 1:M) {
system.time({
    D[[m]] <- purrr::map_dfr(.x = 1:N, .f = ~ calc_dist_par(.x, m=m, locs=locs, grid=grid))
})
# }




# 8.5 secs with 4 workers
# 7.78 secs with 6 workers
# 7.99  secs with 8 workers
D_par <- vector(mode = 'list', length =  M)
# for (m in 1:M) {
# future::plan(multisession, workers = 4)
future::plan(multisession, workers = 6)
# future::plan(multisession, workers = 8)
system.time({
    D_par[[m]] <- furrr::future_map_dfr(.x = 1:N, .f = ~ calc_dist_par(.x, m=m, locs=locs, grid=grid),
                                        .options = furrr_options(seed=NULL))
})
# }


# 19 secs
time_cpp_start <- Sys.time()
D2 <- distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                                   radius = grid$radius[m])
time_cpp_end <- Sys.time()
print(paste("time took:", round(difftime(time_cpp_end, time_cpp_start, units="secs"), digits=2), "seconds"))

all.equal(as.matrix(D_par[[m]][, 1:2]), as.matrix(D2$ind), check.attributes=FALSE)
all.equal(as.matrix(D_par[[m]][, 3]), as.matrix(D2$V), check.attributes=FALSE)
all.equal(as.matrix(D_par[[m]][, 4]), as.matrix(D2$ddistx), check.attributes=FALSE)
all.equal(as.matrix(D_par[[m]][, 5]), as.matrix(D2$ddisty), check.attributes=FALSE)

# this indexing is slower and not currently used
# system.time({
#     D2t <- distance_near_with_ddist_cpp_transpose(t(as.matrix(locs)), t(as.matrix(grid$locs_grid[[m]])),
#                                        radius = grid$radius[m])
# })


#  17.5 sec 5xN as a matrix filled in a loop with byrow=TRUE
system.time({
    D3 <- t(distance_near_loop_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                                   radius = grid$radius[m], byrow=TRUE))
})
# 17.5 sec 5xN as a matrix filled in a loop with byrow=FALSE
system.time({
    D4 <- distance_near_loop_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                                 radius = grid$radius[m], byrow=FALSE)
})

# check equality
all.equal(D3, D4)
all.equal(D2$ind, D3[, 1:2])
all.equal(drop(D2$V), D3[, 3])
all.equal(drop(D2$ddistx), D3[, 4])
all.equal(drop(D2$ddisty), D3[, 5])



# building openmp chunks
# 17 sec for ncores=1
# 17 sec for ncores=2 (running sequentially)
# 17 sec for ncores=4 (running sequentially)
# 17 sec for ncores=5 (running sequentially)
#  sec for ncores=6 (running sequentially)

# # byrow and joint_index play little rowl
# bm <- microbenchmark::microbenchmark(
#     distance_near_chunk_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
#                             radius = grid$radius[m], byrow=TRUE, ncores=6, joint_index=FALSE),
#     distance_near_chunk_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
#                             radius = grid$radius[m], byrow=FALSE, ncores=6, joint_index=FALSE),
#     times = 20
# )
# bm

time_loop_start <- Sys.time()
D6 <- distance_near_loop_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                             radius = grid$radius[m], byrow = FALSE)
time_loop_end <- Sys.time()


# nchunks doesn't seem to matter too much... at least for N approx 65K.
time_openmp_start <- Sys.time()
D5 <- distance_near_chunk_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                              radius = grid$radius[m], byrow=FALSE,
                              nchunks = 6,
                              ncores=6, joint_index=TRUE)
time_openmp_end <- Sys.time()
print(paste("openmp time took:", round(difftime(time_openmp_end, time_openmp_start, units="secs"), digits=2), "seconds"))
print(paste("loop time took:", round(difftime(time_loop_end, time_loop_start, units="secs"), digits=2), "seconds"))
print(paste("cpp time took:", round(difftime(time_cpp_end, time_cpp_start, units="secs"), digits=2), "seconds"))

# check the output
D_out <- do.call(rbind, D5)
all.equal(D2$ind, D_out[, 1:2])
all.equal(drop(D2$V), D_out[, 3])
all.equal(drop(D2$ddistx), D_out[, 4])
all.equal(drop(D2$ddisty), D_out[, 5])


fun <- function(n, nchunks=n) {
    tmp=distance_near_chunk_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                            radius = grid$radius[m], byrow=FALSE, nchunks=nchunks, ncores=n, joint_index=TRUE)
    rm(tmp)
}
fun0 <- function() {
    tmp=distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(grid$locs_grid[[m]]),
                             radius = grid$radius[m])
    rm(tmp)
}

bm <- microbenchmark::microbenchmark(
    fun0(),
    # fun(1),
    # fun(2),
    fun(4),
    fun(6),
    fun(6, nchunks=100*6),
    fun(8),
    fun(10),
    fun(16),
    fun(24),
    fun(36),
    fun(48),
    # fun(56),
    times=10)
bm
print(bm, unit='relative')
autoplot(bm)


# when ncores = 1
# all.equal(D2$ind, D5[[1]][, 1:2])
# all.equal(drop(D2$V), D5[[1]][, 3])
# all.equal(drop(D2$ddistx), D5[[1]][, 4])
# all.equal(drop(D2$ddisty), D5[[1]][, 5])




# microbenchmark::microbenchmark()



