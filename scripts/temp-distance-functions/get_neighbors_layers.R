#' Calculate the neighbors from a coarse gridcell to a finer gridcell
#'
#' @param grid_cell The index (coarser grid cell id)
#' @param g The graph between the coarser and finer grid
#'
#' @return A tibble containing two columns, the gridcell grid_cell and the neighbor grid_cell_neighbors
#' @export
#'
#' @examples
get_neighbors_layers <- function(grid_cell, g) {
    grid_cell_neighbors <- neighbors(g, grid_cell, mode = "out")
    return(tibble(grid_cell_neighbors = as.integer(grid_cell_neighbors), grid_cell = grid_cell))
}
