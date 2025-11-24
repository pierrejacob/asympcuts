#' Function for performing Bayesian linear regression (Metropolis-within-Gibbs);
#' there are two parameters: beta (vector) and sigma 
#' @param final_N final number of iterations (after burn-in and thinning)
#' @param y Outcome variable
#' @param x Matrix with covariates
#' @param prior_sd standard deviation on the normal prior
#' @param sigma_prop standard deviation on the proposal on sigma
#' @param thin use every thin-th iteration
#' @param burn_in number of burn-in iterations (same for both modules)
#' @param initial_state inintial state of the chain; if NULL then the MLE is used

linear_regression_mcmc <- function(N_final, 
                                   x,
                                   y,
                                   prior_sd,
                                   sigma_prop, 
                                   thin =1,
                                   burn_in =1000,
                                   scale = 1){
  N <- N_final*thin + burn_in
  p <- ncol(x)
  n <- nrow(x)
  freq_model <- lm(y~x-1) 
  beta <- freq_model %>% coefficients
  sigma <- sqrt(freq_model$residuals %>% var())
  Sigma_prop <- scale*(sigma^2)*solve(t(x)%*%x)
  
  log_target <- function(beta, sigma){
    log_lik <- dnorm(y - as.vector(x%*%beta),
                     mean = 0, sd = sigma, log = TRUE) %>% sum()
    prior_beta <- sum(dnorm(beta, mean = 0, sd = prior_sd, log = TRUE))
    prior_sigma <- dgamma(sigma, 0.01, 0.01, log = TRUE)
    log_lik + prior_beta + prior_sigma
  }
  log_pdf <- log_target(beta, sigma)
  
  accept_beta <- rep(FALSE, N)
  accept_sigma <- rep(FALSE, N)
  
  chain_beta <- matrix(NA, ncol = p, nrow = N)
  chain_sigma <- rep(NA, N)
  
  for(i in 1:N){
    proposed_beta <- MASS::mvrnorm(1, beta, Sigma_prop)
    proposed_pdf <- log_target(proposed_beta, sigma)
    
    log_u <- log(runif(1))
     if(log_u < (proposed_pdf - log_pdf)){
       accept_beta[i] <- TRUE
       log_pdf <- proposed_pdf
       beta <- proposed_beta
     }
    
    proposed_sigma <- rnorm(1, sigma, sigma_prop)
    proposed_pdf <- log_target(beta, proposed_sigma)
    
    if(log_u < (proposed_pdf - log_pdf)){
      accept_sigma[i] <- TRUE
      log_pdf <- proposed_pdf
      sigma <- proposed_sigma
    }
    
    chain_beta[i,] <- beta
    chain_sigma[i] <- sigma
  }
  
  return_vec <- seq(from = burn_in, by = thin, length.out = N_final)
  chain_sigma <- chain_sigma[return_vec]
  chain_beta <- chain_beta[return_vec,]
  return(list(chain_beta = chain_beta,
         chain_sigma = chain_sigma,
         accept_beta = accept_beta,
         accept_sigma = accept_sigma))
}