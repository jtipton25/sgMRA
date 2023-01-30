calc_r_grid <- function(grid) {

    M <- length(grid$locs_grid)
    r_grid  <- rep(0, M)

    for (m in 1:M) {
        # assumes a rectangular grid
        delta_x <- grid$delta_x[[m]]
        delta_y <- grid$delta_y[[m]]
        r_grid[m] <- sqrt(delta_x^2 + delta_y^2)
    }

    return(r_grid)
}
