#' Generate data.
#'
#' @param n Number of trials.
#' @param n_time Number of time points.
#' @param theta Named list of parameter values.
generate_data <- function(n, n_time, theta) {
  x <- seq(0, 1, length.out = n_time)
  K_f <- sq_exp_kernel(x, theta$rho, nugget = 1e-6)
  K_f_inv <- solve(K_f)

  f <- MASS::mvrnorm(n = 1, rep(0, n_time), K_f)
  y <- matrix(NA, nrow = n_time, ncol = n)
  z <- matrix(NA, nrow = n_time, ncol = n)
  mu <- matrix(NA, nrow = n_time, ncol = n)
  for (i in 1:n) {
    K_i <- get_K_i(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
    mu[, i] <- K_i %*% K_f_inv %*% f
    z[, i] <- arima.sim(model = list(ar = theta$phi), sd = theta$sigma, n = n_time)
    y[, i] <- mu[, i] + z[, i]
  }
  return(list(y = y, f = f, z = z, mu = mu))
}

#' Get Cov(y_i, f)
#'
#' @param x sequence of observation times
#' @param theta list with named entries rho and tau
get_K_i <- function(x, theta) {
  n_time <- length(x)
  K <- matrix(NA, n_time, n_time)
  for (i in 1:n_time) {
    for (j in 1:n_time) {
      K[i, j] <- exp(-theta$rho^2 / 2 * (x[i] - x[j] - theta$tau)^2)
    }
  }
  return(theta$beta * K)
}

#' Get y hat
#'
#' @param i Subject index (scalar, 1 < i < n).
#' @param f Vector of values for f.
#' @param theta Named list of parameter values.
#' @param K_f_inv Inverse covariance matrix of f.
get_y_hat <- function(i, f, theta, K_f_inv) {
  x <- seq(0, 1, length.out = length(f))
  K_i <- get_K_i(x, list(tau = theta$tau[i], beta = theta$beta[i], rho = theta$rho))
  mu <- K_i %*% K_f_inv %*% f
  return(mu)
}

#' Get y_hat for all trials, output in matrix form
#'
#' @param y Matrix of observed trial data.
#' @param f Vector of f values.
#' @param theta Parameter values.
#' @param K_f_inv Inverse covariance matrix of f (n_time x n_time).
get_y_hat_matrix <- function(y, f, theta, K_f_inv) {
  y_hat <- matrix(nrow = nrow(y), ncol = ncol(y))

  for (i in 1:ncol(y_hat)) {
    y_hat[, i] <- get_y_hat(i, f, theta, K_f_inv)
  }
  return(y_hat)
}

#' Generate covariance matrix for square exponential kernel.
#'
#' @param x Vector of time points.
#' @param rho Length scale.
#' @param alpha Amplitude.
#' @param nugget Covariance nugget.
sq_exp_kernel <- function(x, rho, alpha = 1, nugget = 0.0) {
  K <- toeplitz(alpha^2 * exp(-rho^2 / 2 * x^2))
  diag(K) <- diag(K) + nugget
  return(K)
}


#' Get covariance AR process (Section 4.6.2)
#' @param phi autoregressive parms (p x 1)
#' @param sigma error variance
#' @param n_time number of time points
get_Sigma_nu <- function(phi, sigma, n_time)
{
  Sigma_nu_U <- matrix(0, n_time, n_time)
  tmp <- tacvfARMA(phi = phi, maxLag = n_time, sigma2 = sigma^2)
  for (tt in 1:n_time) {
    Sigma_nu_U[tt, tt:(n_time)] = tmp[1:(n_time-tt+1)]
  }
  Sigma_nu <- ultosymmetric(Sigma_nu_U)
  return(Sigma_nu)
}

#' Get Sigma_y_i (Section 4.6.2)
#' @param beta_i Trial specific amplitude beta_i
#' @param K_f Covariance matrix of f. (n_time x n_time)
#' @param Sigma_nu Covariance AR process (n_time x n_time)
get_Sigma_y_i <- function(beta_i, K_f, Sigma_nu)
{
  return((beta_i^2)*K_f + Sigma_nu)
}

#' Get Sigma_y_i f (Section 4.6.2)
#' @param i trial specific amplitude beta_i
#' @param x Kernel GP
#' @param theta Named list of parameter values
#' @param K_f Covariance matrix of f. (n_time x n_time).
#' @param K_f_inv Inverse covariance matrix of f (n_time x n_time).
#' @param Sigma_nu Covariance AR process (n_time x n_time).
getSigma_y_i_f <- function(i, x, theta, K_f, K_f_inv, Sigma_nu)
{
  Sigma_y_i <- get_Sigma_y_i(theta$beta[i], K_f, Sigma_nu)
  K_i <- get_K_i(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
  Sigma_i <- Sigma_y_i - t(K_i)%*%K_f_inv%*%K_i
  Sigma_i <- (Sigma_i + t(Sigma_i))/2
  return(Sigma_i)
}

#' Summary output MCMC
#' @param results output from fit_RPAGP
#' @param dat_trials simulated data in long format 
#' @param y data
#' @param burn_in burn in 
getSummaryOutput <- function(results, dat_trials, y, burn_in)
{
  n <- dim(y)[2]
  n_time <- dim(y)[1]
  n_iter <- length(results$chain)
  n_final <- length((burn_in+1):n_iter)
  
  y_hat <- getSingleTrialEstimates(results, burn_in)
  probs = c(0.025, 0.5, 0.975)
  
  y_hat_quantiles <- array(NA, dim= c(3, n_time, n))
  for (ii in 1:n) {
    y_hat_quantiles[, , ii] = sapply(1:n_time, function(t) quantile(y_hat[t, ii, ], 
                                                                    probs = probs))
  }
  lower <- reshape2::melt(y_hat_quantiles[1, , ], varnames = c("time", "trial"))$value
  median <- reshape2::melt(y_hat_quantiles[2, , ], varnames = c("time", "trial"))$value
  upper <- reshape2::melt(y_hat_quantiles[3, , ], varnames = c("time", "trial"))$value
  
  out <- dat_trials %>% mutate(lwr = lower, med = median, upr = upper)
  
  return(tibble(out))
}

#' get singleTrialEstimates 
#' @param results output from fit_RPAGP
#' @param burn_in burn_in period
getSingleTrialEstimates <- function(results, burn_in)
{
  n_iter <- length(results$chain)
  n_final <- length((burn_in+1):n_iter)
  
  # - beta 
  chain_beta_burned = matrix(NA, n, n_final)
  ss <- 1
  for (tt in (burn_in+1):n_iter) {
    chain_beta_burned[, ss] <- results$chain[[tt]]$beta 
    ss <- ss + 1
  }
  # - f 
  chain_f_burned = matrix(NA, n_time, n_final)
  ss <- 1
  for (tt in (burn_in+1):n_iter) {
    chain_f_burned[, ss] <- results$chain_f[[tt]]
    ss <- ss + 1
  }
  # - 
  y_hat = array(NA, dim = c(n_time, n, n_final))
  for (tt in 1:n_final) {
    for (ii in 1:n) {
      y_hat[, ii, tt] = chain_beta_burned[ii, tt] * chain_f_burned[, tt]
    }
  }
  return(y_hat)
}


#' Transform an upper triangular matrix to symmetric
#' @param m upper triangular matrix
ultosymmetric <- function(m) 
{
  m = m + t(m) - diag(diag(m))
  return (m)
}

