#' Function for performing Bayesian logistic regression
#' @param final_N final number of iterations (after burn-in and thinning)
#' @param y Outcome variable
#' @param x Matrix with covariates
#' @param log_prior function taking beta and evaluating log prior
#' @param beta_scale the proposal covariance matrix is given by the approx.
#' cov matrix at the MLE times 2.38^2/dimension; but this can be further multipled
#' by beta_scale
#' @param burn_in number of burn-in iterations (same for both modules)
#' @param thin use every thin-th iteration
#' @param initial_state inintial state of the chain; if NULL then the MLE is used

logistic_regression_mcmc <- function(final_N,
                                     y,
                                     x,
                                     log_prior,
                                     beta_scale = 1,
                                     burn_in = 1000,
                                     thin = 1,
                                     initial_state = NULL){
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(n == length(y))

  # Find MLE
  mle_point <- (glm(y ~ x -1, family="binomial")$coefficients) %>% unname

  # Find an approximate covariance matrix of the posterior
  score <- (x%*%mle_point)  %>% as.vector
  w_vector <- exp(score)/(1+ exp(score))^2
  w_matrix <- diag(w_vector)

  # Use optimal scaling formula
  proposal_matrix <- beta_scale*2.38^2*solve(t(x)%*%w_matrix%*%x)/p

  # Define the total number of iterations
  n_iter <- burn_in + thin*final_N

  # Draw vectors that will be used later for proposal
  proposed_vectors <- mvtnorm::rmvnorm(n_iter, rep(0, p), proposal_matrix)

  # Set the start values and objects to store results
  chain <- matrix(NA, nrow = n_iter, ncol = p)
  accept <- rep(FALSE, n_iter)
  accept_prob <- rep(FALSE, n_iter)
  if(is.null(initial_state)){
    initial_state <- mle_point + mvtnorm:: rmvnorm(1,
                                                   rep(0, p),
                                                   proposal_matrix)[1,]
  }
  w1 <- rep(1, n)
  current_state <- initial_state
  current_density <- log_prior(current_state) + weighted_logistic_log_likC(current_state,
                                                                     y, x, w1, p, n)
  # Run the main MCMC loop
  for(i in 1:n_iter){
    # Propose a new state
    proposed_state <- current_state + proposed_vectors[i,]

    # Calculating proposed density
    proposed_density <- log_prior(proposed_state) + weighted_logistic_log_likC(proposed_state,
                                                                         y, x, w1, p, n)
    # Accept/reject step
    log_acc_prob <- proposed_density - current_density

    if(log(runif(1)) < log_acc_prob){
      current_state <- proposed_state
      current_density <- proposed_density
      accept[i] <- TRUE
    }
    chain[i,] <- current_state

  }

  result <- list(chain = chain[seq(from = burn_in+1, by = thin, to = n_iter),],
                 accept = accept,
                 proposal_matrix = proposal_matrix)
  return(result)
}

