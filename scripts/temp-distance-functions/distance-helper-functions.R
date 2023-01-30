# helper functions for hierarchical distances





# expand_matrix_list2 <- function(A_list) {
#     M <- length(A_list)
#     n <- rep(0, M)
#     p <- rep(0, M)
#
#
#     rows <- vector(mode = 'list', length = M)
#
#     for (m in 1:M) {
#         n[m] <- nrow(A_list[[m]])
#         p[m] <- ncol(A_list[[m]])
#     }
#
#     for (m in 1:M) {
#
#         if (m == 1) {
#             C <-  Matrix(0, nrow = n[m], ncol = sum(p[2:M]), sparse = TRUE)
#             rows[[m]] <- cbind(A_list[[m]], C)
#
#         } else if (m == M) {
#             B <-  Matrix(0, nrow = n[m], ncol = sum(p[1:(M-1)]), sparse = TRUE)
#             rows[[m]] <- cbind(B, A_list[[m]])
#         } else {
#             B <-  Matrix(0, nrow = n[m], ncol = sum(p[1:(m-1)]), sparse = TRUE)
#             C <-  Matrix(0, nrow = n[m], ncol = sum(p[(m+1):M]), sparse = TRUE)
#             rows[[m]] <- cbind(B, cbind(A_list[[m]], C))
#         }
#     }
#     all_rows <- do.call(rbind, rows)
#     A <- expand_matrix(all_rows)
#     return(A)
# }






