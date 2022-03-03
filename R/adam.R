#' ADAM stochastic process update
#'
#' @param k The gradient descent iteration
#' @param grad A list of gradient vectors
#' @param m A list of momentums
#' @param v A list of velocities
#' @param beta1
#' @param beta2
#' @param epsilon
#'
#' @return
#' @export
#'
#'
adam <- function(k, grad, m, v, beta1=0.9, beta2=0.999, epsilon=1e-8) {
    # iteration k
    # learning rate (eta)
    # current gradient grad
    # current momentum m
    # current velocity v
    # beta1
    # beta2
    # epsilon

    m_hat <- m
    v_hat <- v
    # update m and v
    for (j in 1:length(m)) {
        m[[j]] <- beta1 * m[[j]] + (1 - beta1) * grad[[j]]
        v[[j]] <- beta2 * v[[j]] + (1 - beta2) * grad[[j]]^2
        m_hat[[j]] <- m[[j]] / (1 - beta1^k)
        v_hat[[j]] <- v[[j]] / (1 - beta2^k)
    }

    # update m_hat and v_hat


    return(list(m_hat = m_hat, v_hat = v_hat, m = m, v = v, epsilon = epsilon))
}
