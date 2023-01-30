
get_nearest_gridcell_from_graph <- function(i, locs, grid, g) {
    return(tibble(i = i, idx = trace_neighbors_from_graph(locs[i, ], grid, g)))
}
