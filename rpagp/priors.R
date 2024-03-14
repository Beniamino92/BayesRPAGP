prior <- list()

#' Prior for tau (latency).
#'
#' @param tau tau value.
#' @param tau_prior_sd SD of prior distribution.
prior[["tau"]] <- function(tau, tau_prior_sd) {
  result <- 0
  n <- length(tau)
  Sigma <- tau_prior_sd^2 * (diag(1, n) - matrix(1 / (n+1), n, n))
  result <- result + mvtnorm::dmvnorm(tau[1:(n)], rep(0, n),
                                      Sigma, log = T)
  return(result)
}


#' Prior for rho (GP length scale).
#'
#' @param rho rho value.
#' @param rho_prior_shape Shape parameter of prior.
#' @param rho_prior_scale Scale parameter of prior.
prior[["rho"]] <- function(rho, rho_prior_shape, rho_prior_scale) {
  dgamma(rho, shape = rho_prior_shape, scale = rho_prior_scale, log = T)
}
