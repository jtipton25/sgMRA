get_pairwise_distance <- function(i, g, nearest_idx, locs, trunc) {
    dat_neighbors <- get_neighbors(i, g, nearest_idx)
    if (nrow(dat_neighbors) <= 1) {
        return(tibble(i = numeric(), j = numeric(), D = numeric()))
    } else {
        pairs <- as_tibble(t(Rfast::comb_n(dat_neighbors$v_idx, 2))) %>%
            rename(i = V1, j = V2)
        if (nrow(pairs) == 1) {
            pairs$D <- sqrt(sum((locs[pairs$i, ] - locs[pairs$j, ])^2))
        } else {
            pairs$D <- sqrt(rowSums((locs[pairs$i, ] - locs[pairs$j, ])^2))
        }
        # drop the points below the threshold
        pairs <- pairs %>%
            filter(D <= trunc)
        return(pairs)
    }
}
