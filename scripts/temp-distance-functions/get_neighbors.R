#' Get closest gridcell neighbor
#'
#' Calculate the pairwise neighbors given the absolute closest gridpoint to each location in nearest_idx
#'
#' @param grid_cell The index (coarser grid cell id)
#' @param g The graph between the coarser and finer grid
#' @param nearest_idx # A vector of length equal to the number of observation points that represents the nearest gridcell to the observation location
#'
#' @return A tibble containing two columns, the gridcell grid_cell and the pairwise neighbor grid_cell_neighbors
#' @export
#'
#' @examples
get_neighbors <- function(grid_cell, g, nearest_idx) {
    grid_cell_neighbors <- which(nearest_idx %in% neighbors(g, grid_cell, mode = "in"))
    return(tibble(grid_cell = grid_cell,  grid_cell_neighbors = grid_cell_neighbors))
    # should this really be the following?
    # return(tibble(grid_cell = grid_cell,  grid_cell_neighbors = grid_cell_neighbors))
}
