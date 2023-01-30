euclid_distance <- function(locs, locs_grid) {
    sqrt(rowSums((locs - locs_grid)^2))
}


