dwendland_basis <- function(d, radius) {
    d_rad <- d / radius
    1 / radius * (- 56 / 3 * d_rad * (5 * d_rad + 1) * (1 - d_rad)^5) * (d_rad < 1)
}
