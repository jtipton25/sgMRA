

expand_matrix_symmetric <- function(A) {
    m <- nrow(A)
    n <- ncol(A)
    B <- Matrix(0, nrow = m, ncol = m, sparse = TRUE)
    C <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
    cbind(rbind(B, t(A)), rbind(A, C))
    # rbind(A, Matrix(0, nrow = n-m, ncol = n, sparse = TRUE))
}
