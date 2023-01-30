
trace_neighbors_from_graph <- function(locs, grid, g) {

    if ((length(locs) != 2) | !is.matrix(locs))
        stop("locs must be a matrix with 1 row and 2 columns")

    # coarse_idx <- which.min(sqrt((locs[1] - grid$locs_grid[[1]][, 1])^2 +
    #                                  (locs[2] - grid$locs_grid[[1]][, 2])^2))
    coarse_idx <- which.min(fields::rdist(locs, grid$locs_grid[[1]]))
    idx <- neighbors(g, coarse_idx, mode = 'all')
    return(idx)
}
