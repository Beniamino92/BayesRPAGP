#' Fit RPAGP Model
#'
#' @param y data
#' @param n_iter number of posterior samples
#' @param theta0 parameter initializations
#' @param hyperparam hyperparameters
#' @param pinned_point pinned point for identifiability of structural signal
#' @param pinned_value value of f at the pinned point
#'
fit_rpagp <- function(y, n_iter, theta0, hyperparam, pinned_point, pinned_value = 1) 
{
  chain <- vector(mode = "list", length = n_iter)
  chain_f <- vector(mode = "list", length = n_iter)
  chain_y_hat <- vector(mode = "list", length = n_iter)
  chain_z <- vector(mode = "list", length = n_iter)
  x <- seq(0, 1, length.out = nrow(y))
  z <- matrix(0, nrow = n_time, ncol = n)
  p <- length(theta0$phi)
  
  chain[[1]] <- theta0
  chain_f[[1]] <- sample_f(y, chain[[1]], n_draws = 1, nugget = 1e-6)
  K_f <- sq_exp_kernel(x, chain[[1]]$rho, nugget = 1e-9)
  K_f_inv <- solve(K_f)
  start <- Sys.time()
  
  
  for (iter in 2:n_iter) {
    
    if((iter %% (n_iter/10)) == 0) {
      cat(" ...", as.integer((iter/n_iter)*100), "% \n")
    }
    
    current <- chain[[iter-1]]
    
    # Sample f and rescale
    f  <- sample_f(y, current, n_draws = 1)
    f <- pinned_value * f / f[pinned_point]
    
    # Update y hat
    y_hat <- get_y_hat_matrix(y, f, current, K_f_inv)
    
    # Sample betas
    current$beta <- sample_beta(y = y,
                                current = current,
                                y_hat = y_hat,
                                hyperparam = hyperparam)
    
    # Sample taus and rho
    f <- sample_f(y, current, n_draws = 1)
    f <- pinned_value * f / f[pinned_point]
    current$tau <- sample_tau(y, f, current, hyperparam, K_f, K_f_inv)
    current$rho <- sample_rho(y, f, current, hyperparam)
    K_f <- sq_exp_kernel(x, current$rho, nugget = 1e-6)
    K_f_inv <- solve(K_f)
    
    # Update y hat
    y_hat <- get_y_hat_matrix(y, f, current, K_f_inv)
    
    # Compute residuals and sample residual parameters
    z <- y - y_hat
    ar_post <- sample_AR(z, p)
    
    # Record draws from current iteration
    chain_f[[iter]] <- f
    chain[[iter]] <- current 
    chain[[iter]]$phi <- as.vector(ar_post$rho)
    chain[[iter]]$sigma <- sqrt(ar_post$sigma2)
    chain_y_hat[[iter]] <- y_hat
    chain_z[[iter]] <- z
  }
  cat("\n")
  end <- Sys.time()
  print(end - start)
  
  return(list(chain = chain,
              chain_f = chain_f,
              chain_y_hat = chain_y_hat,
              chain_z = chain_z, 
              runtime = (end-start)))
}
