
#' @note Posterior Bootstrap for both stages of the epidemiological example;
#' @note at the first stage we allow for using Posterior Bootstrap with
#' @note pseudosamples or with prior penalization (with w_0=1); note that in
#' @note any case we don't need to run optimization because in this case we
#' @note have the exact formula for the argmax
#' @param N number of samples
#' @param v_0 value of v_0
#' @param nhpv vector with numbers of HPV cases (as above)
#' @param Npart vector with numbers representing poulations (as above)
#' @param ncases vector with numbers of cases (as above)
#' @param Npop_normalized vector with numbers of follow up years (as above)
#' @param starting_point starting point for optimization for the second stage
#' @param pseudosamples TRUE if the pseudosamples method should be used at the first stage,
#' FALSE otherwise
#' @param control list of parameters to be passed to the optimizer
#' @param nr_cores number of cores

epidemiological_full_pb <- function(N,
                                    v_0,
                                    nhpv,
                                    Npart,
                                    ncases,
                                    Npop_normalized,
                                    starting_point,
                                    pseudosamples = FALSE,
                                    control = list(reltol = 1e-10,
                                                   maxit = 1000),
                                    nr_cores = 6){
  n2 <- length(ncases)
  result_list <- mclapply(1:N, function(i){
    if (pseudosamples) {
      theta_1 <- sapply(1:n2, function(k){
        prob_values <- rbeta(2, 1, 1)
        pseudo_success <- rbinom(1, 1, prob_values[1]) + rbinom(1, 1, prob_values[2])
        weights_success <- sum(rexp(nhpv[k] + pseudo_success, 1))
        weights_failure <- sum(rexp(Npart[k] - nhpv[k] + 2- pseudo_success, 1))
        weights_success/(weights_failure+weights_success)
      })
    } else {
      theta_1 <- sapply(1:n2, function(k){
        weights <- rexp(Npart[k], 1)
        w2 <- sum(weights[(nhpv[k]+1):Npart[k]])
        w1 <- sum(weights) - w2
        (w1+1)/(w1+w2+2)
      })
    }

    weights <- rexp(n2, 1)
    function_to_optim <- function(theta_2){
      weighted_log_lik <- epidemiology_weighted_loglikC(theta_1,
                                                        theta_2,
                                                        weights,
                                                        ncases,
                                                        Npop_normalized)
      log_prior  <- dnorm(theta_2, 0, sqrt(1000), log = TRUE) %>% sum()
      -weighted_log_lik - v_0*log_prior
    }

    current_optim <- tryCatch(optim(par =  starting_point,
                                    fn = function_to_optim,
                                    method = 'BFGS',
                                    control = control),
                              error = function(e)  return(list(error = TRUE)))
    if (isTRUE(current_optim$error)){
      theta_2 <- rep(NA, 2)
    } else {
      if (current_optim$convergence == 0){
        theta_2 <- current_optim$par
      } else {
        theta_2 <- rep(NA, 2)
      }
    }

    list(theta_1 = theta_1, theta_2 = theta_2)
  }, mc.cores = nr_cores)

  theta_1_df <- lapply(result_list, function(x) x[[1]]) %>% unlist %>%
    matrix(ncol = n2, byrow = TRUE) %>% as.data.frame()
  theta_2_df <- lapply(result_list, function(x) x[[2]]) %>% unlist %>%
    matrix(ncol = 2, byrow = TRUE) %>% as.data.frame()

  list(theta_1_df = theta_1_df,
       theta_2_df = theta_2_df)
}


N <- 10000
v_0 <- 1
starting_point <- c(-1,14)

# Posterior Bootstrap with prior penalization at both stages
result_pb_penal <- epidemiological_full_pb(N,
                                           v_0,
                                           nhpv,
                                           Npart,
                                           ncases,
                                           Npop_normalized,
                                           starting_point,
                                           pseudosamples = FALSE,
                                           control = list(reltol = 1e-10,
                                                          maxit = 1000))

pbmi_theta1 <- result_pb_penal$theta_1_df
pbmi_theta2 <- result_pb_penal$theta_2_df

save(pbmi_theta1, pbmi_theta2, file = "inst/epidemio/epidemio_pbmi.RData")
