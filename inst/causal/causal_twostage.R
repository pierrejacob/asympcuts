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

## import functions for MCMC and PBMI
source('inst/causal/causal_inference_linear_regression_mcmc.R')
source('inst/causal/causal_inference_logistic_regression_mcmc.R')
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

#### Cut Bayesian approach ####
N1 <- 1000
stage_1_mcmc <- logistic_regression_mcmc(final_N = N1,
                                         y = lalonde$treat,
                                         x = stage_1_model_mat,
                                         log_prior = normal_log_prior,
                                         beta_scale =1,
                                         burn_in = 5000,
                                         thin = 20)
nr_quantiles <-5
res_mcmc_list <- mclapply(1:N1, function(i){
  beta_1 <- stage_1_mcmc$chain[i,]
  prop_scores_factor <- as.vector(stage_1_model_mat%*%beta_1) %>% gtools::quantcut(q = nr_quantiles)
  x_part_2 <- model.matrix(~prop_scores_factor) %>%
    as.data.frame() %>%
    as.matrix %>%
    unname()
  mat_part_2 <- cbind(lalonde$treat, x_part_2)
  colnames(mat_part_2) <- c('treatment', paste0('V', 1:nr_quantiles))

  stage2_mcmc <- linear_regression_mcmc(N_final = 10,
                                        x = mat_part_2 ,
                                        y = standardize(lalonde$re78),
                                        prior_sd = 1000,
                                        sigma_prop = 0.01,
                                        thin =200,
                                        burn_in =3000,
                                        scale = 1)
  return(stage2_mcmc$chain_beta)

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

#### Cut-Laplace approximation ####

laplace_cut_results <- laplace_cut(logposterior1, logposterior2, colMeans(stage_1_mcmc$chain), colMeans(df_mcmc[,1:6]))
laplace_cut_results$asympvar[11,11]

save(laplace_cut_results, file = 'inst/causal/causal_twostep_laplace_results.RData')

ggplot(df_mcmc, aes(x = V1, y = ..density..)) + geom_density() +
  geom_density(data = data.frame(V1 = rnorm(1e5, laplace_cut_results$thetahat[11], sqrt(laplace_cut_results$asympvar[11,11]))))

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
                                                  nr_cores = 6,
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
                                                        nr_cores = 6,
                                                        control = list(reltol = 1e-10, maxit = 1000),
                                                        same_weights = FALSE)
df_pb <- pb_res$beta_2_df %>% as.data.frame()
df_pb$method <- 'PBMI'
colnames(df_pb)

df_pb_refresh <- pb_res_refresh$beta_2_df %>% as.data.frame()
df_pb_refresh$method <- 'PB refresh'

df <- rbind(df_mcmc, df_pb, df_pb_refresh)
#df
saveRDS(df, file ='inst/causal/causal_twostage_results.rds')
#
# # Visualization to compare PBMI and cut Bayesian
# ggplot(df, aes(x=V1, col = method)) + geom_density()
# ggplot(df, aes(x=V2, col = method)) + geom_density()
# ggplot(df, aes(x=V3, col = method)) + geom_density()
# ggplot(df, aes(x=V4, col = method)) + geom_density()
# ggplot(df, aes(x=V5, col = method)) + geom_density()
# ggplot(df, aes(x=V6, col = method)) + geom_density()



