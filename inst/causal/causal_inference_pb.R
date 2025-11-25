library(Rcpp)
library(RcppArmadillo)
library(optimx)
library(gtools)
library(parallel)
library(dplyr)

#'@param N number of samples
#'@param x matrix of covariates (including intercept)
#'@param treatment the treatment variable
#'@param outcome the outcome variable
#'@param prior_sds_1 standard deviations of priors for module 1
#'@param prior_sds_2 standard deviations of priors for module 2
#'@param ps_continuous TRUE if consider propensity scores as a continuous variable
#'and FALSE otherwise; if FALSE we cut into quantiles
#' @param nr_quantiles if ps_continuous FALSE, how many quantiles are used
#' @param w_0 non-negative value calibrating the impact of the prior on theta_1
#' @param v_0 non-negative value calibrating the impact of the prior on theta_2
#' @param nr_cores number of cores used
#' @param control list of arguments to be passed to the optimizer
#' @param same_weights TRUE if we use Algorithm 2, FALSE if we use ALgorithm 1
#'

cut_model_logistic_linear_reg_bootstrap <- function(N,
                                                    x,
                                                    treatment,
                                                    outcome,
                                                    prior_sds_1,
                                                    prior_sds_2,
                                                    ps_continuous = FALSE,
                                                    nr_quantiles = 5,
                                                    w_0 = 1,
                                                    v_0 = 1,
                                                    nr_cores = 1,
                                                    control = list(reltol = 1e-10, maxit = 1000),
                                                    same_weights = TRUE){
  n <- nrow(x)
  p <- ncol(x)

  result_list <-  mclapply(1:N, function(iter){
    # ------------ First part of the cut model
    weights <- rexp(n, 1)
    function_to_optim <- function(beta){
      weighted_log_lik <- weighted_logistic_log_likC(beta, treatment, x, weights, p, n)
      log_prior <- sum(sapply(1:p, function(j) dnorm(beta[j],
                                                     mean = 0,
                                                     sd = prior_sds_1[j],
                                                     log = TRUE)))
      -weighted_log_lik - w_0*log_prior
    }
    start_beta <- rep(0, p)
    optim_part_1 <- tryCatch(optim(par = start_beta,
                                   fn = function_to_optim,
                                   method = 'BFGS',
                                   control = control),
                             error = function(e)  return(list(error = TRUE)))
    if(isTRUE(optim_part_1$error)){
      result <- list(beta_1 = rep(NA, p),
                     beta_2 = rep(NA, 1 + nr_quantiles))
      return(result)
    } else {
      if(optim_part_1$convergence == 0){
        beta_1 <- optim_part_1$par
      } else {
        result <- list(beta_1 = rep(NA, p),
                       beta_2 = rep(NA, 1 + nr_quantiles))
        return(result)
      }
    }
    # ------------ Second part of the cut model
    # Define second model matrix based on propensity scores
    if(ps_continuous){
      x_part_2 <- cbind(as.vector(x%*%beta_1), rep(1, n))
    }else{
      prop_scores_factor <- as.vector(x%*%beta_1) %>% gtools::quantcut(q = nr_quantiles)
      x_part_2 <- model.matrix(~prop_scores_factor) %>%
        as.data.frame() %>%
        as.matrix %>%
        unname()
    }
    x_part_2 <- cbind(treatment, x_part_2)
    p_2 <- ncol(x_part_2)

    # defining the new target function
    if(isFALSE(same_weights)){
      weights <- rexp(n, 1)
    }

    function_to_optim <- function(theta){
      beta <- theta[1:p_2]
      sigma <- theta[p_2+1]
      if (sigma <= 0) return(-Inf)
      weighted_log_lik <- weighted_linear_regressionC(beta, outcome, x_part_2, weights, p_2, n)
      log_prior <- sum(sapply(1:p_2, function(j) dnorm(beta[j],
                                                       mean = 0,
                                                       sd = prior_sds_2[j],
                                                       log = TRUE)))
      -weighted_log_lik - v_0*(log_prior + dgamma(sigma, 0.01, 0.01, log = TRUE))
    }
    start_beta <- rep(0, p_2)
    optim_part_2 <- tryCatch(optim(par = c(start_beta, 1),
                                   fn = function_to_optim,
                                   method = 'BFGS',
                                   control = control),
                             error = function(e)  return(list(error = TRUE)))

    if(isTRUE(optim_part_2$error)){
      theta_2 <- rep(NA, p_2+1)

    } else {
      if(optim_part_2$convergence == 0){
        theta_2 <- optim_part_2$par
      } else {
        theta_2 <- rep(NA, p_2+1)
      }
    }
    result <- list(beta_1 = beta_1,
                   theta_2 = theta_2)
    return(result)

  }, mc.cores = nr_cores)

  beta_1_df <- do.call('rbind', lapply(result_list, function(x) x$beta_1))
  colnames(beta_1_df) <- paste0('V', 1:ncol(beta_1_df))
  theta_2_df <- do.call('rbind', lapply(result_list, function(x) x$theta_2))
  colnames(theta_2_df) <- paste0('V', 1:ncol(theta_2_df))

  list(beta_1_df = beta_1_df,
       theta_2_df = theta_2_df)
}




