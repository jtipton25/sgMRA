# Binary space partitioning https://en.wikipedia.org/wiki/Binary_space_partitioning
# KD tree https://en.wikipedia.org/wiki/K-d_tree

library(Matrix)
library(fields)
library(tidyverse)

set.seed(101)
# create a 1-d grid
N_coarse = 2^8
N_fine = 2^10


# create an lattice graph
# g = Graphs.SimpleGraphs.grid([10, 4])
# draw(PNG("wheel10.png", 16cm, 16cm), gplot(g))

grid_coarse = seq(-1, 1, length.out=N_coarse)
grid_fine = seq(-1, 1, length.out=N_fine)

# observations
N = 50000
obs = runif(N, -1, 1)

D_obs = rdist(obs)
D_obs_coarse = rdist(obs, grid_coarse)
D_coarse_fine = rdist(grid_coarse, grid_fine)
D_fine = rdist(grid_fine)

radius = 0.0005
# Construct a graph to link coarse grid to fine grid
# nearest coarse grid point to each fine point
nearest_obs_to_coarse = apply(D_obs_coarse, 1, which.min)
nearest_coarse_to_fine = apply(D_coarse_fine, 1, which.min)
nearest_neighbors_fine = apply(D_fine, 2, function(x) which(x <= radius))


find_nearest_fine = function(obs, nearest_obs_to_coarse, nearest_coarse_to_fine, nearest_neighbors_fine,
                             D_obs_coarse, D_coarse_fine, D_fine) {
    # top-down tree
    # start with observation, end with fine grid neighbors
    # Question: how do I make this tree bottom up?
    closest_coarse = rep(NA, N)
    closest_fine = rep(NA, N)
    obs_idx = list()
    fine_neighbors = list()
    for (i in 1:N) {
        closest_coarse[i] = nearest_obs_to_coarse[i]
        closest_fine[i] = nearest_coarse_to_fine[closest_coarse[i]]
        fine_neighbors[[i]] = unlist(nearest_neighbors_fine[closest_fine[i]])
        obs_idx[[i]] = rep(i, length(fine_neighbors[[i]]))
    }
    neighbors_dict = list(
        "obs_idx" = obs_idx,
        "closest_coarse" = closest_coarse,
        "closest_fine" = closest_fine,
        "fine_neighbors" = fine_neighbors
    )
    return (neighbors_dict)
}

neighbors_dict = find_nearest_fine(
    obs,
    nearest_obs_to_coarse,
    nearest_coarse_to_fine,
    nearest_neighbors_fine,
    D_obs_coarse,
    D_coarse_fine,
    D_fine
)

obs_vec = c(unlist(neighbors_dict[["obs_idx"]]))
neighbors_vec = c(unlist(neighbors_dict[["fine_neighbors"]]))
neighbors_mat = sparseMatrix(i=obs_vec, j=neighbors_vec, x=rep(1, length(neighbors_vec)))

# pairwise distance over graph
for (i in 1:N) {
    idx1 = which(neighbors_mat[i, ] == 1)
    idx_neighbors = neighbors_mat[, idx1, drop=FALSE]@i
}


# for observation i
plot_neighbors <- function(i, make_plot=TRUE) {
    idx1 = which(neighbors_mat[i, ] == 1)
    idx_neighbors = neighbors_mat[, idx1, drop=FALSE]@i+1

    dat <- data.frame(obs = obs) %>%
        mutate(neighbor = case_when(D_obs[i, ] <= radius ~ "neighbor",
                                    TRUE ~ "not neighbor"))
    dat_grid = data.frame(grid_fine = grid_fine)

    message("All neighbors ",
            ifelse(all(which(dat$neighbor == "neighbor") %in% idx_neighbors), "are", "are not"),
            " in the hierarchical subset")

    if (make_plot) {
        ggplot() +
            geom_point(data=dat, aes(x = obs, y = 0.0), alpha = 0.2) +
            geom_point(data=dat_grid, aes(x = grid_fine, y = 0.0),
                       color = "red", alpha = .2) +
            geom_point(data = NULL, aes(obs[idx_neighbors], y=0.), color = "green") +
            geom_point(data = NULL, aes(obs[i], y=0.), color = "blue") +
            geom_point(data = dat, aes(x = obs, y=0.05, color = neighbor), alpha = 0.1)
    }
}

plot_neighbors(sample(N, 1), make_plot = TRUE)


# for i in 1:N_fine
# print(i)
# end

# Now I have the neighboring fine scale observations, I need to find all
# points that share the closest fine neighbor with the neighbor index





# neighborhood for coarse grid
g_coarse = Graphs.SimpleGraphs.grid([N_coarse, N_coarse])
g_fine = Graphs.SimpleGraphs.grid([N_fine, N_fine])
# draw(PNG("g_coarse.png", 16cm, 16cm), gplot(g_coarse))
# draw(PNG("g_fine.png", 16cm, 16cm), gplot(g_fine))

# for i in 1:(N_coarse+N_fine)
#     add_vertex!.(g)
# end

# for i in 1:N_coarse
#     add_edge!.(g, i, nearest_idx[i])
# end




