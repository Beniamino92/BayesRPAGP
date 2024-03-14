
#' Likelihood
#'
#' @param y Matrix of observed trial data (n_time x n).
#' @param f Vector of values for f (n_time).
#' @param theta Named list of parameter values.
#' @param K_f Matrix, Kernel GP
#' @param K_f_inv Matrix, inv Kernel GP
likelihood <- function(y, f, theta, K_f, K_f_inv) {
  n <- ncol(y)
  n_time <- nrow(y)
  p <- length(theta$phi)
  
  x <- seq(0, 1, length.out = n_time)
  result <- 0
  z <- matrix(0, nrow = n_time, ncol = n)
  eps <- matrix(0, nrow = n_time, ncol = n)
  Sigma_nu <- get_Sigma_nu(theta$phi, theta$sigma, n_time)
  
  tmp <- 0.0
  for (i in 1:n) {
    Sigma_y_i_f <- getSigma_y_i_f(i, x, theta, K_f, K_f_inv, Sigma_nu)
    Sigma_y_i_f <- (Sigma_y_i_f + t(Sigma_y_i_f))/2
    mu <- get_y_hat(i, f, theta, K_f_inv = K_f_inv)
    tmp <- tmp +  dmvnorm(y[, i] , mean = mu, sigma = Sigma_y_i_f, log = TRUE)
  }
  
  out = ifelse(tmp == -Inf, -1e10, tmp)
  return(out)
}
