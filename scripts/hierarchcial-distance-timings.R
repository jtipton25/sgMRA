# references:
# https://perso.esiee.fr/~aigerd/PID2967435.pdf
# https://www.sciencedirect.com/science/article/pii/S1877705817341747
# https://dl.acm.org/doi/fullHtml/10.1145/3488377
#

# naive attempt
# Goal: given a set of locations and an MRA grid, find all points that are within a given radius of each of the finest gridpoints and calculate these pairwise distances

library(patchwork)
library(tidyverse)
library(Matrix)
library(igraph)
library(fastmatch)
library(data.table)

# load the experimental R functions
file_list <- list.files(here::here("scripts", "temp-distance-functions"), full.names = TRUE)
invisible(lapply(file_list, source))
Rcpp::sourceCpp("~/sgMRA/src/dist_near_cpp.cpp")


# timings of functions ----

# naive attempt
# Goal: given a set of locations and an MRA grid, find all points that are within a given radius of each of the finest gridpoints and calculate these pairwise distances

library(patchwork)
library(tidyverse)
library(Matrix)
library(igraph)
library(fastmatch)

# load the experimental R functions
file_list <- list.files(here::here("scripts", "temp-distance-functions"), full.names = TRUE)
invisible(lapply(file_list, source))


set.seed(1)
n <- c(50, 100, 500, 1000, 5000, 10000)
dat_timings <- data.frame(N = n,
                          time_setup = rep(0, length(n)),
                          time_graph = rep(0, length(n)),
                          time_direct = rep(0, length(n)))
for (N in n) {
    message("Fitting sample size ", N)
    M <- 3
    trunc <- 0.01
    n_coarse_grid <- 10
    locs <- matrix(runif(N*2), N, 2)

    start_setup <- Sys.time()
    grid <- sgMRA::make_grid(locs,
                             min_x = min(locs[, 1]), max_x = max(locs[, 1]),
                             min_y = min(locs[, 2]), max_y = max(locs[, 2]),
                             M = M, n_coarse_grid = n_coarse_grid, n_padding = 0L)

    # build the graphs for each resolution -- pretty sure this is not needed
    r_grid <- calc_r_grid(grid)



    # for each observation, calculate the nearest fine grid cell
    D <- fields::rdist(locs, grid$locs_grid[[M]])
    fine_idx <- apply(D, 1, which.min)

    # create a graph of neighbors among the finest resolution of the grid
    D_grid <- fields::rdist(grid$locs_grid[[M]])
    A <- Matrix((D_grid < 1.01 * r_grid[M]), sparse = TRUE) * 1
    A <- expand_matrix(A)
    g <- igraph::graph_from_adjacency_matrix(A)
    neighbors = ego(g, order = 1, nodes = fine_idx)

    end_setup <- Sys.time()
    # loc_neighbors <- vector(mode = 'list')

    # Need to make this MUCH MUCH MUCH faster -- can I use the graph traversal instead?
    # This still seems to be O(N^2)
    # for (i in 1:N) {
    #     loc_neighbors[[i]] <- which(sapply(neighbors, function(x) fine_idx[i] %fin% x))
    # }
    # calc_pairwise_distance <- function(i, locs, loc_neighbors) {
    #     j <- loc_neighbors[[i]]
    #     if (length(j) == 0) {
    #         return(tibble(i = numeric(0), j = numeric(0), D = numeric(0)))
    #     }
    #     return(tibble(i = i, j = j, D = drop(fields::rdist(matrix(locs[i, ], 1, 2), matrix(locs[j, ], length(j), 2)))))
    # }
    #
    # # Need to make this MUCH faster
    # tmp <- vector(mode = 'list', length = N)
    # for (i in 1:N) {
    #     tmp[[i]] <- calc_pairwise_distance(i, locs, loc_neighbors)
    # }
    # # Need to make this MUCH faster
    # dat_distance <- map_dfr(.x = 1:N,
    #                         .f = calc_pairwise_distance,
    #                         locs = locs,
    #                         loc_neighbors = loc_neighbors)
    #

    # turn the neighbors into a data.frame

    start_graph <- Sys.time()

    n_neighbors <- unlist(lapply(neighbors, function(x) length(x)))
    dat_all <- tibble(i = rep(1:N, n_neighbors), j = unlist(neighbors))
    # i=111
    # all.equal(dat_all$i[dat_all$j == fine_idx[i]], which(sapply(neighbors, function(x) fine_idx[i] %fin% x)))

    calc_pairwise_distance2 <- function(i, locs, dat_all) {
        j <- dat_all$i[dat_all$j == fine_idx[i]]
        if (length(j) == 0) {
            return(tibble(i = numeric(0), j = numeric(0), D = numeric(0)))
        }
        D <- c(fields::rdist(matrix(locs[i, ], 1, 2), matrix(locs[j, ], length(j), 2)))
        zero_idx <- which(D == 0)
        if (length(zero_idx) > 0) {
            D <- D[-zero_idx]
            j <- j[-zero_idx]
        }
        return(tibble(i = i, j = j, D = D))
    }

    # Need to make this MUCH faster
    # tmp <- vector(mode = 'list', length = N)
    # for (i in 1:N) {
    #     tmp[[i]] <- calc_pairwise_distance2(i, locs, dat_all)
    # }
    # Need to make this MUCH faster
    dat_distance <- map_dfr(.x = 1:N,
                            .f = calc_pairwise_distance2,
                            locs = locs,
                            dat_all = dat_all)
    dat_distance <- dat_distance[dat_distance$D <= trunc, ]
    D_locs <- sparseMatrix(i = dat_distance$i, j = dat_distance$j, x = dat_distance$D, dims = c(N, N))

    end_graph <- Sys.time()

    # Need to figure out why these don't agree
    start_direct <- Sys.time()
    D_locs_direct <- fields::rdist(locs)
    D_locs_direct[D_locs_direct > trunc] <- 0
    D_locs_sparse <- Matrix(D_locs_direct, sparse = TRUE)
    end_direct <- Sys.time()

    start_cpp <- Sys.time()
    D_cpp <- distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(locs),
                                          radius = trunc, n_neighbors = 200)
    D_cpp_sparse <- sparseMatrix(i = D_cpp$ind[, 1], j = D_cpp$ind[, 2], x = c(D_cpp$V), dims = c(N, N))
    end_cpp <- Sys.time()



    dat_timings$time_setup[dat_timings$N == N] <- difftime(end_setup, start_setup, units = "secs")
    dat_timings$time_graph[dat_timings$N == N] <- difftime(end_graph, start_graph, units = "secs")
    dat_timings$time_direct[dat_timings$N == N] <- difftime(end_direct, start_direct, units = "secs")
    dat_timings$time_cpp[dat_timings$N == N] <- difftime(end_cpp, start_cpp, units = "secs")

}


dat_plot <- dat_timings %>%
    pivot_longer(cols = c("time_setup", "time_graph", "time_direct", "time_cpp"),
                 names_to = "operation",
                 values_to = "time")
ggplot(dat_plot, aes(x = N, y = time, color = operation)) +
    geom_line()


# proof of concept with visualization ----

set.seed(1)
M <- 3
N <- 5000
trunc <- 0.01
n_coarse_grid <- 10
locs <- matrix(runif(N*2), N, 2)
grid <- sgMRA::make_grid(locs,
                         min_x = min(locs[, 1]), max_x = max(locs[, 1]),
                         min_y = min(locs[, 2]), max_y = max(locs[, 2]),
                         M = M, n_coarse_grid = n_coarse_grid, n_padding = 0L)

# build the graphs for each resolution -- pretty sure this is not needed
r_grid <- calc_r_grid(grid)



# for each observation, calculate the nearest fine grid cell
D <- fields::rdist(locs, grid$locs_grid[[M]])
fine_idx <- apply(D, 1, which.min)

# create a graph of neighbors among the finest resolution of the grid
D_grid <- fields::rdist(grid$locs_grid[[M]])
A <- Matrix((D_grid < 1.01 * r_grid[M]), sparse = TRUE) * 1
A <- expand_matrix(A)
g <- igraph::graph_from_adjacency_matrix(A)
neighbors = ego(g, order = 1, nodes = fine_idx)
# loc_neighbors <- vector(mode = 'list')

# Need to make this MUCH MUCH MUCH faster -- can I use the graph traversal instead?
# This still seems to be O(N^2)
# for (i in 1:N) {
#     loc_neighbors[[i]] <- which(sapply(neighbors, function(x) fine_idx[i] %fin% x))
# }
# calc_pairwise_distance <- function(i, locs, loc_neighbors) {
#     j <- loc_neighbors[[i]]
#     if (length(j) == 0) {
#         return(tibble(i = numeric(0), j = numeric(0), D = numeric(0)))
#     }
#     return(tibble(i = i, j = j, D = drop(fields::rdist(matrix(locs[i, ], 1, 2), matrix(locs[j, ], length(j), 2)))))
# }
#
# # Need to make this MUCH faster
# tmp <- vector(mode = 'list', length = N)
# for (i in 1:N) {
#     tmp[[i]] <- calc_pairwise_distance(i, locs, loc_neighbors)
# }
# # Need to make this MUCH faster
# dat_distance <- map_dfr(.x = 1:N,
#                         .f = calc_pairwise_distance,
#                         locs = locs,
#                         loc_neighbors = loc_neighbors)
#

# turn the neighbors into a data.frame
n_neighbors <- unlist(lapply(neighbors, function(x) length(x)))


# tibble dat_all
system.time({
    dat_all <- tibble(i = rep(1:N, n_neighbors), j = unlist(neighbors))

    for (i in 1:5000) {
        j <- dat_all$i[dat_all$j == fine_idx[i]]
    }
})

# data.frame dat_all
system.time({
    dat_all <- data.frame(i = rep(1:N, n_neighbors), j = unlist(neighbors))

    for (i in 1:5000) {
        j <- dat_all$i[dat_all$j == fine_idx[i]]
    }
})


# data.table dat_all
system.time({
    dat_all <- data.table(i = rep(1:N, n_neighbors), j = unlist(neighbors))

    for (i in 1:5000) {
        j <- dat_all$i[dat_all$j == fine_idx[i]]
    }
})

# i=111
# all.equal(dat_all$i[dat_all$j == fine_idx[i]], which(sapply(neighbors, function(x) fine_idx[i] %fin% x)))

calc_pairwise_distance2 <- function(i, locs, dat_all, type = "data.frame") {
    # data.frame and data.table are fastest. tibble is slowest
    j <- dat_all$i[dat_all$j == fine_idx[i]]
    if (length(j) == 0) {
        if (type == "tibble") {
            return(tibble(i = numeric(0), j = numeric(0), D = numeric(0)))
        } else if (type == "data.frame") {
            return(data.frame(i = numeric(0), j = numeric(0), D = numeric(0)))
        } else if (type == "data.table") {
            return(data.table(i = numeric(0), j = numeric(0), D = numeric(0)))
        }
    }
    D <- c(fields::rdist(matrix(locs[i, ], 1, 2), matrix(locs[j, ], length(j), 2)))
    zero_idx <- which(D == 0)
    if (length(zero_idx) > 0) {
        D <- D[-zero_idx]
        j <- j[-zero_idx]
    }
    if (type == "tibble") {
        return(tibble(i = i, j = j, D = D))
    } else if (type == "data.frame") {
        return(data.frame(i = i, j = j, D = D))
    } else if (type == "data.table") {
        return(data.table(i = i, j = j, D = D))
    }
}

# Need to make this MUCH faster
# tmp <- vector(mode = 'list', length = N)
# for (i in 1:N) {
#     tmp[[i]] <- calc_pairwise_distance2(i, locs, dat_all)
# }
# Need to make this MUCH faster
dat_distance <- map_dfr(.x = 1:N,
                        .f = calc_pairwise_distance2,
                        locs = locs,
                        dat_all = dat_all)
dat_distance <- dat_distance[dat_distance$D <= trunc, ]
D_locs <- sparseMatrix(i = dat_distance$i, j = dat_distance$j, x = dat_distance$D, dims = c(N, N))


# Need to figure out why these don't agree
D_locs_direct <- fields::rdist(locs)
D_locs_direct[D_locs_direct > trunc] <- 0
D_locs_sparse <- Matrix(D_locs_direct, sparse = TRUE)

D_cpp <- distance_near_with_ddist_cpp(as.matrix(locs), as.matrix(locs),
                                      radius = trunc, n_neighbors = 200)
D_cpp_sparse <- sparseMatrix(i = D_cpp$ind[, 1], j = D_cpp$ind[, 2], x = c(D_cpp$V), dims = c(N, N))

str(D_locs)
str(D_locs_sparse)
str(D_cpp_sparse)

D_locs[dat_distance$i[1:5], dat_distance$j[1:5]]
D_locs_sparse[dat_distance$i[1:5], dat_distance$j[1:5]]
D_cpp_sparse[dat_distance$i[1:5], dat_distance$j[1:5]]



i <- 340
dat_locs <- data.frame(x = locs[, 1], y = locs[, 2], grid_cell = 1:N) %>%
    # mutate(color = ifelse(grid_cell == i, "location", ifelse(grid_cell %in% loc_neighbors[[i]], "neigbor", "exterior")))
    mutate(color = ifelse(grid_cell == i, "location", ifelse(grid_cell %in% dat_all$i[dat_all$j == fine_idx[i]], "neigbor", "exterior")))

ggplot(dat_locs) +
    geom_point(aes(x, y, color = color), size = 0.2) +
    geom_point(data = grid$locs_grid[[M]], aes(x = Var1, y = Var2), alpha = 0.2, size = 0.1) +
    scale_color_viridis_d() +
    ggtitle("Hierarchical search area")


# calculate the hierarchical graph -- use each layer to successively build the others
# generate a "graph" between the first two grid layers ----

A <- vector(mode = 'list', length = M-1)
g <- vector(mode = 'list', length = M-1)
nearest_idx <- vector(mode = 'list', length = M-1)
neighbors <- vector(mode = 'list', length = M-1)
dat_neighbors <- vector(mode = 'list', length = M-1)


D <- fields::rdist(grid$locs_grid[[1]], grid$locs_grid[[2]])
# Adjacency matrix between layers of the grid
A[[1]] <- Matrix((D < 1.01 * r_grid[1]), sparse = TRUE) * 1
A[[1]] <- expand_matrix(A[[1]])
g[[1]] <- igraph::graph_from_adjacency_matrix(A[[1]])
nearest_idx[[1]] <- apply(D, 2, which.min)
dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[2]]),
                              .f = get_neighbors,
                              g = g[[1]],
                              nearest_idx = nearest_idx[[1]])

# dat_n_layres=map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
#                                   get_neighbors_layers,
#                                   g = g[[1]]) %>%
#         distinct() %>%
#         mutate(layer = paste(1, 2, sep = '~'))

if (M > 2) {
    for (m in 2:(M-1)) {

        D <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
        A[[m]] <- Matrix((D <= 1.01 * r_grid[m]), sparse = TRUE) * 1
        A[[m]] <- expand_matrix(A[[m]])
        # adjacency graph between grid layer m and grid layer m+1
        g[[m]] <- igraph::graph_from_adjacency_matrix(A[[m]])

        nearest_idx[[m]] <- apply(D, 2, which.min)
        # get the neighbors between layers
        dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
                                      .f = get_neighbors,
                                      g = g[[m]],
                                      nearest_idx = nearest_idx[[m]])
    }
}

expect_snapshot(dat_neighbors)

# visualize the neighbors at each layer
# layer 2 neighbors for a layer 1 gridcell
i=44
dat_grid <- as_tibble(grid$locs_grid[[1]]) |>
    rename(x = Var1, y = Var2) |>
    mutate(grid_cell = 1:nrow(grid$locs_grid[[1]])) %>%
    mutate(color = ifelse(grid_cell == i, "center", "exterior"))

grid_fine_idx <- dat_neighbors[[1]] %>%
    filter(grid_cell == i) %>%
    pull(grid_cell_neighbors)

dat_grid_fine <- as_tibble(grid$locs_grid[[2]]) |>
    rename(x = Var1, y = Var2) |>
    mutate(grid_cell = 1:nrow(grid$locs_grid[[2]])) %>%
    mutate(color = ifelse(grid_cell %fin% grid_fine_idx, "neighbor", "fine grid")) |>
    mutate(alpha = ifelse(grid_cell %fin% grid_fine_idx, 0.75, 0.25))

p_grid <- ggplot() +
    geom_point(data = dat_grid_fine, aes(x, y, color = color, alpha = alpha)) +
    geom_point(data = dat_grid, aes(x, y, color = color), inherit.aes = TRUE) +
    scale_color_viridis_d(begin = 1, end = 0)
p_grid


# layer 3 neighbors for a layer 2 gridcell
i=78
dat_grid <- as_tibble(grid$locs_grid[[2]]) |>
    rename(x = Var1, y = Var2) |>
    mutate(grid_cell = 1:nrow(grid$locs_grid[[2]])) %>%
    mutate(color = ifelse(grid_cell == i, "center", "exterior"))

grid_fine_idx <- dat_neighbors[[2]] %>%
    filter(grid_cell == i) %>%
    pull(grid_cell_neighbors)
dat_grid_fine <- as_tibble(grid$locs_grid[[3]]) |>
    rename(x = Var1, y = Var2) |>
    mutate(grid_cell = 1:nrow(grid$locs_grid[[3]])) %>%
    mutate(color = ifelse(grid_cell %fin% grid_fine_idx, "neighbor", "fine grid")) |>
    mutate(alpha = ifelse(grid_cell %fin% grid_fine_idx, 0.75, 0.25))

p_grid <- ggplot() +
    geom_point(data = dat_grid_fine, aes(x, y, color = color, alpha = alpha)) +
    geom_point(data = dat_grid, aes(x, y, color = color), inherit.aes = TRUE) +
    scale_color_viridis_d(begin = 1, end = 0)
p_grid
