

get_nearest_gridcell <- function(i, locs, grid, dat_neighbors,  use_tidy = FALSE, use_fast_match = TRUE) {
    return(tibble(i = i, idx = trace_neighbors(locs=matrix(locs[i, ], 1, 2), grid, dat_neighbors, use_tidy, use_fast_match)))
}
