## load packages
library(asympcuts)
library(dplyr)
library(ggplot2)
library(MatchIt)
## load data 'lalonde' from MatchIt
## Description: this is a subsample of the data from the treated group in the
## National Supported Work Demonstration (NSW) and the comparison sample from
## the Population Survey of Income Dynamics (PSID). This data was previously
## analyzed extensively by Lalonde (1986) and Dehejia and Wahba (1999).
data(lalonde)

## compile Rcpp function to evaluate (weighted) log-likelihood of logistic
## regression & linear regression
library(Rcpp)
sourceCpp('inst/causal/weighted_loglik.cpp')

## import functions for PBMI
source('inst/causal/causal_inference_pb.R')

## Standardize the data
## Can we use the 'scale' function in the base package?
standardize <- function(x) {
  (x-mean(x))/sd(x)
}

## have a look
head(lalonde)
lalonde$race <- factor(lalonde$race, levels = c('white', 'hispan', 'black'))
lalonde$age <- standardize(lalonde$age)
lalonde$educ <- standardize(lalonde$educ)
lalonde$re74 <- standardize(lalonde$re74)
lalonde$re75 <- standardize(lalonde$re75)
lalonde$age_sq <- standardize((lalonde$age)^2)

normal_log_prior <- function(theta){
  sum(dnorm(theta, 0, 1000, TRUE))
}

stage_1_model_mat <- model.matrix(treat ~ age + age_sq + educ + race + married + nodegree + re74 + re75,
                                  data = lalonde)

p_1 <- ncol(stage_1_model_mat)
n <- nrow(stage_1_model_mat)
w1 <- rep(1, n)
logposterior1 <- function(theta1) normal_log_prior(theta1) + weighted_logistic_log_likC(theta1, lalonde$treat, stage_1_model_mat, w1, p_1, n)
## Compute Laplace approximation
theta1hat <- (glm(lalonde$treat ~ stage_1_model_mat -1, family="binomial")$coefficients) %>% unname
theta1hat <- optim(par = theta1hat, fn = function(x) -logposterior1(x))$par
print(theta1hat)

## Compute Laplace approximation
asympvar1 <- solve(hessian(func = function(theta) -logposterior1(theta),
                           x = theta1hat))

#### Cut Bayesian approach ####
## using the generic RWMRTH implementation
target <- list(dimension = p_1, density = function(theta1){
  return(sapply(1:dim(theta1)[1], function(ichain) logposterior1(theta1[ichain,])))
})

## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(rep(0.1, p_1))
tuning_parameters$niterations <- 10000
tuning_parameters$nchains <- 4
tuning_parameters$adaptation <- 1000
tuning_parameters$rinit <- function(nchains){
  chains <- matrix(rnorm(nchains*target$dimension, 0, 1), ncol = target$dimension)
  return(chains)
}
## initialization
current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
## MCMC run
rwmrth_res_1 <- rwmrth(target, tuning_parameters)
## default choice of burnin
burnin <- tuning_parameters$niterations/2

mcmc_samples_stage1 <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
  rwmrth_res_1$chains[[ichain]][(burnin+1):tuning_parameters$niterations,]
}

## first stage accuracy
# library(gridExtra)
# laplace1_samples <- fast_rmvnorm(1e6, theta1hat, asympvar1)
#
# df_stage1 <- (rbind(data.frame(method = "mcmc", theta = mcmc_samples_stage1),
#                     data.frame(method = "laplace", theta = laplace1_samples[,1:p_1])))
# g1 <- ggplot(df_stage1, aes(x = theta.1, fill = method)) +
#   geom_density(alpha = 0.5)
# g2 <- ggplot(df_stage1, aes(x = theta.2, fill = method)) +
#   geom_density(alpha = 0.5)
# g3 <- ggplot(df_stage1, aes(x = theta.3, fill = method)) +
#   geom_density(alpha = 0.5)
# grid.arrange(g1, g2, g3,  ncol = 3)

N1 <- 5000
selected_indices <- round(seq(1, nrow(mcmc_samples_stage1), length.out = N1))


nr_quantiles <-5
res_mcmc_list <- mclapply(selected_indices, function(i){
  beta_1 <- mcmc_samples_stage1[i,]
  prop_scores_factor <- as.vector(stage_1_model_mat%*%beta_1) %>% gtools::quantcut(q = nr_quantiles)
  x_part_2 <- model.matrix(~prop_scores_factor) %>%
    as.data.frame() %>%
    as.matrix %>%
    unname()
  mat_part_2 <- cbind(lalonde$treat, x_part_2)
  colnames(mat_part_2) <- c('treatment', paste0('V', 1:nr_quantiles))
  p_2 <- ncol(mat_part_2) + 1
  y = standardize(lalonde$re78)
  x = mat_part_2
  prior_sd = 1000
  freq_model <- lm(y~x-1)
  betahat <- freq_model %>% coefficients
  sigma <- sqrt(freq_model$residuals %>% var())
  Sigma_prop <- (sigma^2)*solve(t(x)%*%x)
  log_target <- function(theta2){
    beta <- theta2[1:(p_2-1)]
    sigma <- theta2[p_2]
    if (sigma <= 0) return(-Inf)
    prior_sigma <- dgamma(sigma, 0.01, 0.01, log = TRUE)
    log_lik <- dnorm(y - as.vector(x %*% matrix(beta, ncol = 1)),
                     mean = 0, sd = sigma, log = TRUE) %>% sum()
    prior_beta <- sum(dnorm(beta, mean = 0, sd = prior_sd, log = TRUE))
    return(log_lik + prior_beta + prior_sigma)
  }

  ## using the generic RWMRTH implementation
  target <- list(dimension = p_2, density = function(theta){
    return(sapply(1:dim(theta)[1], function(ichain) log_target(theta[ichain,])))
  })

  ## tuning parameters for rwmrth
  tuning_parameters <- list()
  cov_proposal <- matrix(0, nrow = p_2, ncol = p_2)
  cov_proposal[1:(p_2-1), 1:(p_2-1)] <- Sigma_prop
  cov_proposal[p_2, p_2] <- 0.1
  tuning_parameters$cov_proposal <- cov_proposal
  tuning_parameters$niterations <- 2000
  tuning_parameters$nchains <- 1
  tuning_parameters$adaptation <- 1000
  tuning_parameters$rinit <- function(nchains){
    chains <- matrix(0, nrow = nchains, ncol = target$dimension)
    chains[,1:(p_2-1)] <- fast_rmvnorm(nchains, betahat, Sigma_prop)
    chains[,p_2] <- rexp(nchains, 1/sigma)
    return(chains)
  }
  ## initialization
  current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
  ## MCMC run
  rwmrth_res_2 <- rwmrth(target, tuning_parameters)
  ## trace plot
  ## default choice of burnin
  # burnin <- tuning_parameters$niterations/2
  # matplot(rwmrth_res_2$chains[[1]], type = 'l')
  mcmc_samples_stage2 <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
    rwmrth_res_2$chains[[ichain]][(tuning_parameters$niterations-1):tuning_parameters$niterations,]
  }
  return(mcmc_samples_stage2)
}, mc.cores = 10)

df_mcmc <- do.call('rbind', res_mcmc_list ) %>% as.data.frame()
df_mcmc$method <- 'cut Bayesian'
dim(df_mcmc)

var(df_mcmc$V1)


logposterior2 <- function(theta1, theta2){
  prop_scores_factor <- as.vector(stage_1_model_mat %*% theta1) %>% gtools::quantcut(q = 5)
  x_part_2 <- model.matrix(~prop_scores_factor) %>%
    as.data.frame() %>%
    as.matrix %>%
    unname()
  mat_part_2 <- cbind(lalonde$treat, x_part_2)
  p_2 <- ncol(mat_part_2)
  prior_part <- sum(sapply(1:p_2, function(j) dnorm(theta2[j], mean = 0, sd = 1000, log = TRUE)))
  weights <- rep(1, n)
  part_module_2 <- weighted_linear_regressionC(theta2, standardize(lalonde$re78), mat_part_2, weights, p_2, n)
  return(prior_part + part_module_2)
}

theta2hat <- optim(par = rep(0, 7), fn = function(x) -logposterior2(theta1hat, x))$par
print(theta2hat)

#### Posterior Bootstrap for Modular Inference ####
pb_res <- cut_model_logistic_linear_reg_bootstrap(N = nrow(df_mcmc),
                                                  x = stage_1_model_mat,
                                                  treatment = lalonde$treat,
                                                  outcome = standardize(lalonde$re78),
                                                  prior_sds_1 = rep(100, ncol(stage_1_model_mat)),
                                                  prior_sds_2 = rep(100, nr_quantiles +1),
                                                  ps_continuous = FALSE,
                                                  nr_quantiles = 5,
                                                  w_0 = 0,
                                                  v_0 = 0,
                                                  nr_cores = 10,
                                                  control = list(reltol = 1e-10, maxit = 1000),
                                                  same_weights = TRUE)

# What would happen if we used independently drawn weights (Algorithm 1)
pb_res_refresh <- cut_model_logistic_linear_reg_bootstrap(N = nrow(df_mcmc),
                                                        x = stage_1_model_mat,
                                                        treatment = lalonde$treat,
                                                        outcome = standardize(lalonde$re78),
                                                        prior_sds_1 = rep(100, ncol(stage_1_model_mat)),
                                                        prior_sds_2 = rep(100, nr_quantiles +1),
                                                        ps_continuous = FALSE,
                                                        nr_quantiles = 5,
                                                        w_0 = 0,
                                                        v_0 = 0,
                                                        nr_cores = 10,
                                                        control = list(reltol = 1e-10, maxit = 1000),
                                                        same_weights = FALSE)
df_pb <- pb_res$theta_2_df %>% as.data.frame()
df_pb$method <- 'PBMI'
colnames(df_pb)

df_pb_refresh <- pb_res_refresh$theta_2_df %>% as.data.frame()
df_pb_refresh$method <- 'PB refresh'

df_twostage <- rbind(df_mcmc, df_pb, df_pb_refresh)
#df
save(theta1hat, theta2hat, df_twostage, file ='inst/causal/causal_twostage_results.RData')
#
# # Visualization to compare PBMI and cut Bayesian
# ggplot(df, aes(x=V1, col = method)) + geom_density()
# ggplot(df, aes(x=V2, col = method)) + geom_density()
# ggplot(df, aes(x=V3, col = method)) + geom_density()
# ggplot(df, aes(x=V4, col = method)) + geom_density()
# ggplot(df, aes(x=V5, col = method)) + geom_density()
# ggplot(df, aes(x=V6, col = method)) + geom_density()



