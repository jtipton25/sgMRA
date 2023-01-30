# hierarchical pairwise distance estimates
# Use recursive grid search to find neighbors between two sets of vectors v1 and v2 ----

library(patchwork)
library(tidyverse)
library(Matrix)
library(igraph)
library(fastmatch)

# load the experimental R functions
file_list <- list.files(here::here("scripts", "temp-distance-functions"), full.names = TRUE)
invisible(lapply(file_list, source))


# Single level, find points nearest to gridcell

# example for the function

N <- 5000
v1 <- matrix(runif(N*2), N, 2)
v2 <- matrix(runif(N*2), N, 2)

trunc = 0.1
n_knots = 20
max_levels = 4


# v1 is a Nx2 vector
# v2 is a Nx2 vector
# trunc is a truncation
range_x <- range(range(v1[, 1]), range(v2[, 1]))
range_y <- range(range(v1[, 2]), range(v2[, 2]))
# loop over the hierarchy
# grid <- vector(mode = 'list', length = max_levels)
# nodes <- vector(mode = 'list', length = max_levels)
# # for(j in 1:max_levels) {
# seq_x <- seq(range_x[1], range_x[2], length.out=n_knots)
# seq_y <- seq(range_y[1], range_y[2], length.out=n_knots)
# delta_x <- abs(seq_x[1] - seq_x[2])
# delta_y <- abs(seq_y[1] - seq_y[2])
# grid[[j]] <- tidyr::expand_grid(seq_x, seq_y) %>%
#     as.matrix()


# start with a single layer, create a recursive layer above

seq_x <- seq(range_x[1], range_x[2], length.out=n_knots)
seq_y <- seq(range_y[1], range_y[2], length.out=n_knots)
delta_x <- abs(seq_x[1] - seq_x[2])
delta_y <- abs(seq_y[1] - seq_y[2])
r_grid <- sqrt(delta_x^2 + delta_y^2)
# g <- igraph::make_lattice(length = n_knots, dim=2)
# A <- igraph::as_adjacency_matrix(g)


grid <- tidyr::expand_grid(seq_x, seq_y) %>%
    as.matrix()
# }
# for (j in 1:max_levels) {
# while(min(delta_x, delta_y) > 2 * trunc) assign the grid cells

# construct adjacency matrix using the grid
D_grid <- fields::rdist(grid)
A_grid <- Matrix((D_grid < 1.01 * r_grid) * 1)
g <- igraph::graph_from_adjacency_matrix(A_grid)


# pairwise distance from vectors to the grid to find nearest gridpoint
D1 <- fields::rdist(v1, grid)
D2 <- fields::rdist(v2, grid)
nearest_idx1 <- apply(D1, 1, which.min)
nearest_idx2 <- apply(D2, 1, which.min)

for (i in 1:n_knots^2) {
    # make this into a furrr function of i?
    # neighbors1 <- Matrix((nearest_idx1 %in% neighbors(g, i)) * 1)
    # neighbors2 <- Matrix((nearest_idx2 %in% neighbors(g, i)) * 1)
    pairwise_to_check1 <- which(nearest_idx1 %in% neighbors(g, i, mode = "all"))
    pairwise_to_check2 <- which(nearest_idx2 %in% neighbors(g, i, mode = "all"))
    dat_check <- expand_grid(x = pairwise_to_check1, y = pairwise_to_check2, grid_cell = i)
}

# explore the check
i=188
# neighbors1 <- Matrix((nearest_idx1 %in% neighbors(g, i)) * 1)
# neighbors2 <- Matrix((nearest_idx2 %in% neighbors(g, i)) * 1)
pairwise_to_check1 <- which(nearest_idx1 %in% neighbors(g, i, mode = "all"))
pairwise_to_check2 <- which(nearest_idx2 %in% neighbors(g, i, mode = "all"))
dat_check <- expand_grid(v1 = pairwise_to_check1, v2 = pairwise_to_check2, grid_cell = i)

D <- fields::rdist(v1, v2)
D[cbind(dat_check$v1, dat_check$v2)]
hist(D[cbind(dat_check$v1, dat_check$v2)])


dat_shading <- dat_check %>%
    pivot_longer(cols = c(v1, v2), names_to = "vector", values_to = "v_idx") %>%
    mutate(shading = TRUE)

dat_grid <- as.data.frame(grid) |>
    rename(x = seq_x, y = seq_y) |>
    mutate(grid_cell = 1:nrow(grid)) %>%
    mutate(color = ifelse(grid_cell == i, "center", ifelse(grid_cell %in% neighbors(g, i, mode = "all"), "neighbor", "exterior")))

p_grid <- ggplot(dat_grid, aes(x, y, color = color)) +
    geom_point() +
    scale_color_viridis_d()

dat_plot <- data.frame(x = c(v1[, 1], v2[, 1]), y = c(v1[, 2], v2[, 2]),
                       v_idx = c(1:nrow(v1), 1:nrow(v2)),
                       vector = factor(c(rep("v1", nrow(v1)), rep("v2", nrow(v2))))) %>%
    left_join(dat_shading, by=c("vector", "v_idx"))

p_points <- ggplot() +
    geom_point(data = dat_plot, aes(x, y, color = vector, alpha = dat_shading), alpha = 0.2) +
    geom_point(data = dat_grid, aes(x = x, y = y, color = color), inherit.aes = FALSE) +
    scale_color_viridis_d()

p_grid / p_points
# proportion of filtered points for pairwise calculation
nrow(dat_check) / N^2 * nrow(grid)






# try and write a general function ----
# TODO: Need to get a hierarchical filtering done
M <- 5
N <- 5000
n_coarse_grid <- 10
locs <- matrix(runif(N*2), N, 2)
grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)
trunc = 0.1


# hierarchical_pairwise_distances <- function(locs, grid, trunc = 0.1) {

    # example for the function
    library(tidyverse)
    library(Matrix)
    library(igraph)
    library(fastmatch)
    library(furrr)



    # start with a single layer, create a recursive layer above

    # QUESTION: How to use the resolution information hierarchcially? How to develop a data structure to
    # exploit this hierarchical structure?
    # M = 1

    # construct adjacency matrix using the grid, can cache these and reuse these
    # Move to eval_basis.R file in the make_basis() function

    # build the graphs for each resolution -- pretty sure this is not needed
    # g <- vector(mode = 'list', length = M)
    r_grid  <- rep(0, M)

    for (m in 1:M) {

    #     # assumes a rectangular grid
        delta_x <- grid$delta_x[[m]]
        delta_y <- grid$delta_y[[m]]
        r_grid[m] <- sqrt(delta_x^2 + delta_y^2)
    #
    #     # grid within a resolution
    #     D_grid <- fields::rdist(grid$locs_grid[[m]])
    #     A_grid <- Matrix((D_grid < 1.01 * r_grid[m]) * 1)
    #     g[[m]] <- igraph::graph_from_adjacency_matrix(A_grid)
    }

    # calculate the hierarchical graph -- use each layer to successively build the others



    # generate a "graph" between the different grid layers ----
    g_layers <- vector(mode = 'list', length = M-1)
    nearest_idx_layers <- vector(mode = 'list', length = M-1)
    neighbors_layers <- vector(mode = 'list', length = M-1)
    dat_neighbors <- vector(mode = 'list', length = M-1)
    # calculate the first layer by "brute force search"
    D_grid_layers <- fields::rdist(grid$locs_grid[[1]], grid$locs_grid[[1+1]])
    A_grid_layers <- Matrix((D_grid_layers < 1.01 * r_grid[1]), sparse = TRUE) * 1
    A_grid_layers <- expand_matrix(A_grid_layers)
    # A_grid_layers <- expand_matrix_symmetric(A_grid_layers)
    # adjacency graph between grid layer 1 and grid layer 2
    g_layers[[1]] <- igraph::graph_from_adjacency_matrix(A_grid_layers)
    # nearest coarse grid cell to each fine grid cell
    nearest_idx_layers[[1]] <- apply(D_grid_layers, 2, which.min)
    # nearest_idx_layers[[1]] <- nrow(grid$locs_grid[[1]]) + apply(D_grid_layers, 2, which.min)
    # neighbors_layers[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1+1]]),
    #                                  .f = get_neighbors,
    #                                  g = g_layers[[1]],
    #                                  nearest_idx = nearest_idx_layers[[1]])

    dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
                                  get_neighbors_layers,
                                  g = g_layers[[1]]) %>%
        distinct() %>%
        mutate(layer = paste(1, 2, sep = '~'))

    if (M > 2) {
        for (m in 2:(M-1)) {

            D_grid_layers <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
            A_grid_layers <- Matrix((D_grid_layers <= 1.01 * r_grid[m]), sparse = TRUE) * 1
            A_grid_layers <- expand_matrix(A_grid_layers)
            # adjacency graph between grid layer 1 and grid layer 2
            g_layers[[m]] <- igraph::graph_from_adjacency_matrix(A_grid_layers)

            # get the neighbors between layers
            dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
                                          get_neighbors_layers,
                                          g = g_layers[[m]]) %>%
                distinct() %>%
                mutate(layer = paste(m, m+1, sep = '~'))
        }
    }


    # create one very large graph rather than many small graphs to compare speed

    A_grid_list <- vector(mode = 'list', length = M-1)
    D_grid_layers <- fields::rdist(grid$locs_grid[[1]], grid$locs_grid[[1+1]])
    A_grid_list[[1]] <- Matrix((D_grid_layers < 1.01 * r_grid[1]), sparse = TRUE) * 1
    # A_grid_layers <- expand_matrix(A_grid_layers)
    # adjacency graph between grid layer 1 and grid layer 2
    # g_layers[[1]] <- igraph::graph_from_adjacency_matrix(A_grid_layers)
    # nearest coarse grid cell to each fine grid cell
    # nearest_idx_layers[[1]] <- apply(D_grid_layers, 2, which.min)
    # nearest_idx_layers[[1]] <- nrow(grid$locs_grid[[1]]) + apply(D_grid_layers, 2, which.min)
    # neighbors_layers[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1+1]]),
    #                                  .f = get_neighbors,
    #                                  g = g_layers[[1]],
    #                                  nearest_idx = nearest_idx_layers[[1]])

    # dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
    #                               get_neighbors_layers,
    #                               g = g_layers[[1]]) %>%
    #     distinct() %>%
    #     mutate(layer = paste(1, 2, sep = '~'))

    if (M > 2) {
        for (m in 2:(M-1)) {

            D_grid_layers <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
            A_grid_list[[m]] <- Matrix((D_grid_layers <= 1.01 * r_grid[m]), sparse = TRUE) * 1
            # A_grid_layers <- expand_matrix(A_grid_layers)
            # adjacency graph between grid layer 1 and grid layer 2
            # g_layers[[m]] <- igraph::graph_from_adjacency_matrix(A_grid_layers)

            # get the neighbors between layers
            # dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
            #                               get_neighbors_layers,
            #                               g = g_layers[[m]]) %>%
            #     distinct() %>%
            #     mutate(layer = paste(m, m+1, sep = '~'))
        }
    }

    # AA <- rapply(A_grid_list, expand_matrix)
    # AA <- expand_matrix(A_grid_list[[1]])




    A_grid_all <- expand_matrix_list(A_grid_list)
    g_all <- igraph::graph_from_adjacency_matrix(A_grid_all)
    neighbors(g_all, 501, mode = "all")
    neighbors(g_all, 1, mode='all')

    # alternative approach -- find nearest coarse grid cell then evaluate all possible neighbors at once

    dat_locs_grid_from_graph <- purrr::map_dfr(.x = 1:nrow(locs),
                                    .f = get_nearest_gridcell_from_graph,
                                    locs = locs,
                                    grid = grid,
                                    g = g_all)


    coarse_idx <- apply(fields::rdist(locs, grid$locs_grid[[1]]), 1, which.min)
    zzz=sapply(1:nrow(locs), function(i) neighbors(g_all, coarse_idx[i], mode = 'all'))

    # dat_neighbors <- map_dfr(.x = 1:nrow(grsid$locs_grid[[2-1]]),
    #                              get_neighbors,
    #                              g = g_layers[[2-1]],
    #                              nearest_idx = nearest_idx_layers[[2-1]])

    # plot the neighbors between layers of the grid ----
    # m <- 2
    # m <- 3
    # m <- 4
    p_grid <- vector(mode = 'list', length = M-1)
    p_points <- vector(mode = 'list', length = M-1)
    for (m in 2:M) { # assumes n_coarse_grid = 10
        if (m == 2) {
            i= 275 # m = 2
        } else if (m == 3) {
            i = 435 # m = 3
        } else if (m == 4) {
            i = 785 # m = 4
        } else  if (m == 5) {
            i = 2485 # m = 5
        }

        dat_neighbors_i <- dat_neighbors[[m-1]] %>%
            filter(grid_cell == i)

        dat_shading <- dat_neighbors_i %>%
            mutate(color = "neighbor")

        dat_grid <- as_tibble(grid$locs_grid[[m-1]]) |>
            rename(x = Var1, y = Var2) |>
            mutate(grid_cell = 1:nrow(grid$locs_grid[[m-1]])) %>%
            mutate(color = ifelse(grid_cell == i, "center", "exterior"))
        # mutate(color = ifelse(grid_cell == i, "center", ifelse(grid_cell %in% neighbors(g_layers[[1]], i, mode = "all"), "neighbor", "exterior")))

        p_grid[[m-1]] <- ggplot(dat_grid, aes(x, y, color = color)) +
            geom_point() +
            scale_color_viridis_d(begin = 0.2, end = 0.8)

        dat_plot <- tibble(x = grid$locs_grid[[m]][, 1], y = grid$locs_grid[[m]][, 2],
                           v_idx = nrow(grid$locs_grid[[m-1]]) + 1:nrow(grid$locs_grid[[m]])) %>%
            left_join(dat_shading, by=c("v_idx")) %>%
            mutate(color = replace_na(color, "inner"))
        # mutate(color = ifelse(grid_cell %in% neighbors(g_layers[[1]], i, mode = "all"), "neighbor", "observation"))
        # mutate(alpha = replace_na(alpha, 0.5), color = "observation")

        p_points[[m-1]] <- ggplot() +
            geom_point(data = dat_grid, aes(x = x, y = y, color = color), alpha = 0.5) +
            geom_point(data = dat_plot, aes(x = x, y = y, color = color), alpha = 0.5) +
            scale_color_viridis_d(begin = 1., end = 0.)
    }
    wrap_plots(p_grid, nrow = 2)
    wrap_plots(p_points, nrow = 2)

    # Trace the graph between layers to a terminal node for each observation location ----




    # this is REALLY SLOW. How can this be sped up?
    dat_locs_grid <- purrr::map_dfr(.x = 1:nrow(locs),
                                    .f = get_nearest_gridcell,
                                    locs = locs,
                                    grid = grid,
                                    dat_neighbors = dat_neighbors)


    # using profvis, at least it appears to be linear in complexity of time...


    tmp <- vector(mode='list', length=nrow(locs))
    for (i in 1:nrow(locs)) {
        tmp[[i]] <- get_nearest_gridcell(i, locs = locs,
                                         grid = grid, dat_neighbors = dat_neighbors)

    }

    DD=fields::rdist(locs, grid$locs_grid[[M]])
    idx_test <- dat_locs_grid %>%
        filter(i == 1) %>%
        pull(idx)
    layout(matrix(1:2, 2, 1))
    hist(DD[1, idx_test])
    hist(DD[1, ])




# }

dat_grid <- as_tibble(grid$locs_grid[[4]]) |>
    rename(x = Var1, y = Var2) |>
    mutate(grid_cell = 1:nrow(grid$locs_grid[[4]])) |>
    mutate(color = ifelse(grid_cell %in% idx, "neighbor", "not neighbor"))

p_grid <- ggplot(dat_grid, aes(x, y, color = color)) +
    geom_point(alpha = 0.5) +
    scale_color_viridis_d(end=0.8) #+
    # geom_point(data = data.frame(x = locs1[1], y = locs1[2]),
    #            aes(x = x, y = y), inherit.aes = FALSE)
p_grid
# }




D_grid_layers <- fields::rdist(grid$locs_grid[[1]], grid$locs_grid[[1+1]])
D_grid_layers[D_grid_layers > r_grid[1]] <- 0

tmp = list()
for (i in 1:nrow(grid$locs_grid[[1+1]])) {
    tmp[[i]] <- get_pairwise_distance(i, g_layers[[1]], nearest_idx_layers[[1]],
                                      grid$locs_grid[[1+1]], r_grid[1])
}

dat_distance <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
                        get_pairwise_distance,
                        g = g_layers[[1]],
                        nearest_idx = nearest_idx_layers[[1]],
                        locs = grid$locs_grid[[1+1]],
                        trunc = r_grid[1])

# remove duplicate rows -- lots of duplicates
dat_distance <- dat_distance %>%
    distinct()

D_sparse_from_layers <- Matrix(D_grid_layers, sparse = TRUE)
D_sparse <- sparseMatrix(i = dat_distance$i,
                         j = dat_distance$j,
                         x = dat_distance$D,
                         dims = c(nrow(grid$locs_grid[[1]]), nrow(grid$locs_grid[[1+1]])))

all.equal(D_sparse_from_layers, D_sparse)


A_grid_layers <- Matrix((D_grid < 1.01 * r_grid[m]) * 1)
g_layers[[m]] <- igraph::graph_from_adjacency_matrix(A_grid_layers)
nearest_idx_layers[[m]] <- apply(D_grid_layers, 2, which.min)
neighbors_layers[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m+1]]),
                                 .f = get_neighbors,
                                 g = g_layers[[m]],
                                 nearest_idx = nearest_idx_layers[[m]])

}




##
## Use the graph to evaluate pairwise distance between the locations
##

# build the graphs for each resolution
g <- vector(mode = 'list', length = M)
r_grid  <- rep(0, M)

for (m in 1:M) {

    # assumes a rectangular grid
    delta_x <- grid$delta_x[[m]]
    delta_y <- grid$delta_y[[m]]
    r_grid[m] <- sqrt(delta_x^2 + delta_y^2)

    # grid within a resolution
    D_grid <- fields::rdist(grid$locs_grid[[m]])
    A_grid <- Matrix((D_grid < 1.01 * r_grid[m]) * 1)
    g[[m]] <- igraph::graph_from_adjacency_matrix(A_grid)
}

# pairwise distance from vectors to the grid to find nearest gridpoint
D <- fields::rdist(locs, grid$locs_grid[[m]])
nearest_idx <- apply(D, 1, which.min)



# can parallelize this if needed...
dat_check <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]), get_neighbors, g = g[[m]], nearest_idx = nearest_idx)
# tmp = list()
# for (i in 1:nrow(grid$locs_grid[[m]])) {
#     tmp[[i]] <- get_pairwise_distance(i, g[[m]], nearest_idx, trunc)
# }
dat_distance <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
                        get_pairwise_distance, g = g[[m]], nearest_idx = nearest_idx,
                        locs = locs, trunc = trunc)

# remove duplicate rows
dat_distance <- dat_distance %>%
    distinct()

# calculate the sparse index
dat_check$D <- sqrt(sum((locs[dat_check$v_idx, ] - grid$locs_grid[[m]][dat_check$grid_cell, ])^2))
dat_check

D_full <- fields::rdist(locs)
D_full[D_full > trunc] <- 0
D_sparse_from_full <- Matrix(D_full, sparse = TRUE)
D_sparse <- sparseMatrix(i = dat_distance$i,
                         j = dat_distance$j,
                         x = dat_distance$D,
                         dims = c(nrow(locs), nrow(locs)),
                         symmetric = TRUE)

all.equal(D_sparse_from_full, D_sparse)

# D_sparse <- sparseMatrix(i = dat_check$v_idx, j = dat_check$grid_cell, x = dat_check$D,
#                         dims = c(nrow(locs), nrow(grid$locs_grid[[m]])))
# explore the check
if (m == 1) {
    i <- 193 # m=1
} else if (m ==2) {
    i <- 525 # m=2
} else if (m == 3) {
    i <- 1274 # m=3
}


dat_check_i <- dat_check %>%
    filter(grid_cell == i)

# neighbors1 <- Matrix((nearest_idx1 %in% neighbors(g, i)) * 1)
# neighbors2 <- Matrix((nearest_idx2 %in% neighbors(g, i)) * 1)
# pairwise_to_check1 <- which(nearest_idx %in% neighbors(g, i, mode = "all"))
# dat_check <- tibble(v_idx = pairwise_to_check1, grid_cell = i)

if (nrow(D) < 1000) {
    D <- fields::rdist(locs, grid$locs_grid[[m]])
    D[cbind(dat_check_i$v_idx, i)]
    layout(matrix(1:2, 2, 1))
    hist(D)
    hist(D[cbind(dat_check_i$v_idx, i)])
}

dat_shading <- dat_check_i %>%
    # pivot_longer(cols = c(v1, v2), names_to = "vector", values_to = "v_idx") %>%
    mutate(alpha = 0.25)

dat_grid <- as_tibble(grid$locs_grid[[m]]) |>
    rename(x = Var1, y = Var2) |>
    mutate(grid_cell = 1:nrow(grid$locs_grid[[m]])) %>%
    mutate(color = ifelse(grid_cell == i, "center", ifelse(grid_cell %in% neighbors(g[[m]], i, mode = "all"), "neighbor", "exterior")))

p_grid <- ggplot(dat_grid, aes(x, y, color = color)) +
    geom_point() +
    scale_color_viridis_d()
p_grid

dat_plot <- tibble(x = locs[, 1], y = locs[, 2], v_idx = 1:nrow(locs)) %>%
    left_join(dat_shading, by=c("v_idx")) %>%
    mutate(alpha = replace_na(alpha, 0.05), color = "observation")

p_points <- ggplot() +
    geom_point(data = dat_plot, aes(x, y, alpha = alpha, color = color)) +
    geom_point(data = dat_grid, aes(x = x, y = y, color = color)) +
    scale_color_viridis_d(begin = 0., end = 0.9) +
    scale_alpha_identity()
p_points

p_grid / p_points
# proportion of filtered points for pairwise calculation
nrow(dat_neighbors[[M-1]]) / N^2 * nrow(grid$locs_grid[[m]])
}
