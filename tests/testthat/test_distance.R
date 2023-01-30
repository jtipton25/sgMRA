# testing distance functions

file_list <- list.files(here::here("scripts", "temp-distance-functions"), full.names = TRUE)
invisible(lapply(file_list, source))

library(tidyverse)
library(fields)
library(Matrix)
library(igraph)
library(fastmatch)


library(testthat)


# test the expand_matrix_list function ----
test_that("expand_matrix_list", {
    A <- Matrix(1:20, 4, 5, sparse = TRUE)
    B <- Matrix(21:100, 4, 10, sparse = TRUE)
    expect_snapshot(expand_matrix_list(list(A, B)))

    A <- Matrix(1:20, 5, 4, sparse = TRUE)
    B <- Matrix(21:100, 4, 10, sparse = TRUE)
    expect_error(expand_matrix_list(list(A, B)), "All matrices in the list must have more columns than rows")
})



# test the expand_matrix function ----
test_that("expand_matrix", {
    A <- Matrix(1:20, 4, 5, sparse = TRUE)
    expect_snapshot(expand_matrix(A))

    # must have more columns than rows
    A <- Matrix(1:20, 5, 4, sparse = TRUE)
    expect_error(expand_matrix(A), "A must have more columns than rows")
})



# test the expand_matrix_symmetric function ----
test_that("expand_matrix_symmetric", {

    A <- Matrix(1:20, 4, 5, sparse = TRUE)
    expect_snapshot(expand_matrix_symmetric(A))
})



# test the trace_neighbors function ----
test_that("trace_neighbors", {

    set.seed(1)
    M <- 3
    N <- 50
    n_coarse_grid <- 10
    locs <- matrix(runif(N*2), N, 2)
    grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)

    # build the graphs for each resolution -- pretty sure this is not needed
    r_grid <- calc_r_grid(grid)

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
    dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
                                  get_neighbors_layers,
                                  g = g[[1]]) %>%
        distinct() %>%
        mutate(layer = paste(1, 2, sep = '~'))

    if (M > 2) {
        for (m in 2:(M-1)) {

            D <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
            A[[m]] <- Matrix((D <= 1.01 * r_grid[m]), sparse = TRUE) * 1
            A[[m]] <- expand_matrix(A[[m]])
            # adjacency graph between grid layer m and grid layer m+1
            g[[m]] <- igraph::graph_from_adjacency_matrix(A[[m]])

            # get the neighbors between layers
            dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
                                          get_neighbors_layers,
                                          g = g[[m]]) %>%
                distinct() %>%
                mutate(layer = paste(m, m+1, sep = '~'))
        }
    }


    expect_snapshot(trace_neighbors(matrix(locs[10, ], 1, 2), grid, dat_neighbors))

})

# test the get_nearest_gridcell function ----
test_that("get_nearest_gridcell", {

    set.seed(1)
    M <- 3
    N <- 50
    n_coarse_grid <- 10
    locs <- matrix(runif(N*2), N, 2)
    grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)

    # build the graphs for each resolution -- pretty sure this is not needed
    r_grid <- calc_r_grid(grid)

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
    dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
                                  get_neighbors_layers,
                                  g = g[[1]]) %>%
        distinct() %>%
        mutate(layer = paste(1, 2, sep = '~'))

    if (M > 2) {
        for (m in 2:(M-1)) {

            D <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
            A[[m]] <- Matrix((D <= 1.01 * r_grid[m]), sparse = TRUE) * 1
            A[[m]] <- expand_matrix(A[[m]])
            # adjacency graph between grid layer m and grid layer m+1
            g[[m]] <- igraph::graph_from_adjacency_matrix(A[[m]])

            # get the neighbors between layers
            dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
                                          get_neighbors_layers,
                                          g = g[[m]]) %>%
                distinct() %>%
                mutate(layer = paste(m, m+1, sep = '~'))
        }
    }




    # A_all <- expand_matrix_list(A)
    # g_all <- igraph::graph_from_adjacency_matrix(A_all)
    # neighbors(g_all, 501, mode = "all")
    # neighbors(g_all, 1, mode='all')

    # alternative approach -- find nearest coarse grid cell then evaluate all possible neighbors at once

    # dat_locs_grid_from_graph <- purrr::map_dfr(.x = 1:nrow(locs),
    #                                            .f = get_nearest_gridcell_from_graph,
    #                                            locs = locs,
    #                                            grid = grid,
    #                                            g = g_all)


    # coarse_idx <- apply(fields::rdist(locs, grid$locs_grid[[1]]), 1, which.min)
    # zzz=sapply(1:nrow(locs), function(i) neighbors(g_all, coarse_idx[i], mode = 'all'))

    # dat_neighbors <- map_dfr(.x = 1:nrow(grsid$locs_grid[[2-1]]),
    #                              get_neighbors,
    #                              g = g_layers[[2-1]],
    #                              nearest_idx = nearest_idx_layers[[2-1]])


    # Trace the graph between layers to a terminal node for each observation location ----




    # this is REALLY SLOW. How can this be sped up?
    dat_locs_grid <- purrr::map_dfr(.x = 1:nrow(locs),
                                    .f = get_nearest_gridcell,
                                    locs = locs,
                                    grid = grid,
                                    dat_neighbors = dat_neighbors)



    dat_locs_grid %>%
        mutate(neighbors = TRUE) %>%
        right_join(
            data.frame(x = grid$locs_grid[[M]][, 1], y = grid$locs_grid[[M]][, 2], idx = 1:nrow(grid$locs_grid[[M]])),
            by = 'idx') %>%
        ggplot(aes(x, y, color = neighbors)) +
        geom_point()
})



#
# # test the get_nearest_gridcell_from_graph function ----
# test_that("get_nearest_gridcell_from_graph", {
#
# set.seed(1)
# M <- 3
# N <- 50
# n_coarse_grid <- 10
# locs <- matrix(runif(N*2), N, 2)
# grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)
#
# # build the graphs for each resolution -- pretty sure this is not needed
# r_grid <- calc_r_grid(grid)
#
# # calculate the hierarchical graph -- use each layer to successively build the others
# # generate a "graph" between the first two grid layers ----
#
# A <- vector(mode = 'list', length = M-1)
# g <- vector(mode = 'list', length = M-1)
# nearest_idx <- vector(mode = 'list', length = M-1)
# neighbors <- vector(mode = 'list', length = M-1)
# dat_neighbors <- vector(mode = 'list', length = M-1)
#
#
# D <- fields::rdist(grid$locs_grid[[1]], grid$locs_grid[[2]])
# # Adjacency matrix between layers of the grid
# A[[1]] <- Matrix((D < 1.01 * r_grid[1]), sparse = TRUE) * 1
# A[[1]] <- expand_matrix(A[[1]])
# g[[1]] <- igraph::graph_from_adjacency_matrix(A[[1]])
# nearest_idx[[1]] <- apply(D, 2, which.min)
# dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
#                               get_neighbors_layers,
#                               g = g[[1]]) %>%
#     distinct() %>%
#     mutate(layer = paste(1, 2, sep = '~'))
#
# if (M > 2) {
#     for (m in 2:(M-1)) {
#
#         D <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
#         A[[m]] <- Matrix((D <= 1.01 * r_grid[m]), sparse = TRUE) * 1
#         A[[m]] <- expand_matrix(A[[m]])
#         # adjacency graph between grid layer m and grid layer m+1
#         g[[m]] <- igraph::graph_from_adjacency_matrix(A[[m]])
#
#         # get the neighbors between layers
#         dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
#                                       get_neighbors_layers,
#                                       g = g[[m]]) %>%
#             distinct() %>%
#             mutate(layer = paste(m, m+1, sep = '~'))
#     }
# }
#
#
#
#
# # A_all <- expand_matrix_list(A)
# # g_all <- igraph::graph_from_adjacency_matrix(A_all)
# # neighbors(g_all, 501, mode = "all")
# # neighbors(g_all, 1, mode='all')
#
# # alternative approach -- find nearest coarse grid cell then evaluate all possible neighbors at once
#
# # dat_locs_grid_from_graph <- purrr::map_dfr(.x = 1:nrow(locs),
# #                                            .f = get_nearest_gridcell_from_graph,
# #                                            locs = locs,
# #                                            grid = grid,
# #                                            g = g_all)
#
#
# # coarse_idx <- apply(fields::rdist(locs, grid$locs_grid[[1]]), 1, which.min)
# # zzz=sapply(1:nrow(locs), function(i) neighbors(g_all, coarse_idx[i], mode = 'all'))
#
# # dat_neighbors <- map_dfr(.x = 1:nrow(grsid$locs_grid[[2-1]]),
# #                              get_neighbors,
# #                              g = g_layers[[2-1]],
# #                              nearest_idx = nearest_idx_layers[[2-1]])
#
#
# # Trace the graph between layers to a terminal node for each observation location ----
#
#
#
#
# # this is REALLY SLOW. How can this be sped up?
# dat_locs_grid <- purrr::map_dfr(.x = 1:nrow(locs),
#                                 .f = get_nearest_gridcell,
#                                 locs = locs,
#                                 grid = grid,
#                                 dat_neighbors = dat_neighbors)
#
#
#
# dat_locs_grid %>%
#     mutate(neighbors = TRUE) %>%
#     right_join(
#         data.frame(x = grid$locs_grid[[M]][, 1], y = grid$locs_grid[[M]][, 2], idx = 1:nrow(grid$locs_grid[[M]])),
#         by = 'idx') %>%
#     ggplot(aes(x, y, color = neighbors)) +
#     geom_point()
#
# })
#
#
# # test the trace_matrix_from_graph function ----
# test_that("trace_matrix_from_graph", {
#
# })
#
#
#
#
# # test the get_pairwise_distance function ----
# test_that("get_pairwise_distance", {
#
# })
#
#
# # test the get_pairwise_layers function ----
# test_that("get_pairwise_layers", {
#
# })
#

# test the get_neighbors function ----
test_that("get_neighbors", {

    # Not quite sure what this function does and why
    set.seed(1)
    M <- 4
    N <- 50
    n_coarse_grid <- 10
    locs <- matrix(runif(N*2), N, 2)
    grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)

    # build the graphs for each resolution -- pretty sure this is not needed
    r_grid <- calc_r_grid(grid)

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

})


# test the get_neighbors_layers function ----
test_that("get_neighbors_layers", {

    # This function finds the neighbors between the coarse grid and the lower level fine grid
    set.seed(1)
    M <- 3
    N <- 50
    n_coarse_grid <- 10
    locs <- matrix(runif(N*2), N, 2)
    grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)

    # build the graphs for each resolution -- pretty sure this is not needed
    r_grid <- calc_r_grid(grid)

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
    dat_neighbors[[1]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[1]]),
                                  get_neighbors_layers,
                                  g = g[[1]]) %>%
        distinct() %>%
        mutate(layer = paste(1, 2, sep = '~'))

    if (M > 2) {
        for (m in 2:(M-1)) {

            D <- fields::rdist(grid$locs_grid[[m]], grid$locs_grid[[m+1]])
            A[[m]] <- Matrix((D <= 1.01 * r_grid[m]), sparse = TRUE) * 1
            A[[m]] <- expand_matrix(A[[m]])
            # adjacency graph between grid layer m and grid layer m+1
            g[[m]] <- igraph::graph_from_adjacency_matrix(A[[m]])

            # get the neighbors between layers
            dat_neighbors[[m]] <- map_dfr(.x = 1:nrow(grid$locs_grid[[m]]),
                                          get_neighbors_layers,
                                          g = g[[m]]) %>%
                distinct() %>%
                mutate(layer = paste(m, m+1, sep = '~'))
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
    i=84
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
})


test_that("calc_r_grid", {
    set.seed(1)
    M <- 2
    N <- 50
    n_coarse_grid <- 10
    locs <- matrix(runif(N*2), N, 2)
    grid <- sgMRA::make_grid(locs, M = M, n_coarse_grid = n_coarse_grid, n_padding = 1L)

    # calcuate the radius for a queens-move neighbor search radius
    expect_snapshot(calc_r_grid(grid))
})
