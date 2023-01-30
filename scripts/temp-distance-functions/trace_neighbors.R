trace_neighbors <- function(locs, grid, dat_neighbors, use_tidy = FALSE, use_fast_match = TRUE) {
    # given a location, trace the neighborhood down to the nearest fine resolution grid cell
    if ((length(locs) != 2) | !is.matrix(locs))
        stop("locs must be a matrix with 1 row and 2 columns")
    # loop over locs
    M <- length(grid$locs_grid)
    D <- fields::rdist(locs, grid$locs_grid[[1]])
    if (M > 1) {
        idx <- which.min(D) # nearest coarse grid cell
    }
    for (m in 2:M) {

        # calculate distance at next layer using the neighborhood
        # tidyverse way is slow
        if (use_tidy) {
            if (use_fast_match) {
                idx <- dat_neighbors[[m-1]] %>%
                    filter(grid_cell %fin% idx) %>%
                    pull(grid_cell_neighbors) - nrow(grid$locs_grid[[m-1]])
            } else {

                idx <- dat_neighbors[[m-1]] %>%
                    filter(grid_cell %in% idx) %>%
                    pull(grid_cell_neighbors) - nrow(grid$locs_grid[[m-1]])
            }
        } else {
            if (use_fast_match) {
                idx <- dat_neighbors[[m-1]]$grid_cell_neighbors[dat_neighbors[[m-1]]$grid_cell %fin% idx]
            } else {
                idx <- dat_neighbors[[m-1]]$grid_cell_neighbors[dat_neighbors[[m-1]]$grid_cell %in% idx]
            }
        }

        # if (m < M) {
            D <- fields::rdist(locs, grid$locs_grid[[m]][idx, ])
            # Don't subset the neighborhood for the final distance
            idx <- idx[which.min(D)]
        # }
    }
    return(idx)
}
