library(tidyverse)
library(asympcuts)
rm(list = ls())
set.seed(1)
theme_set(theme_bw())
theme_update(axis.text.x = element_text(size = 20),
             axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 20, margin=margin(20,0,0,0)),
             axis.title.y = element_text(size = 20, angle = 90, margin = margin(0,20,0,0)),
             panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                             colour = "gray"),
             panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                             colour = "gray"),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15),
             title = element_text(size = 30),
             strip.text = element_text(size = 20),
             strip.background = element_rect(fill="white"),
             legend.position = "bottom")


standardize <- function(x) {
  (x-mean(x))/sd(x)
}
## compile Rcpp function to evaluate (weighted) log-likelihood of logistic
## regression & linear regression
library(Rcpp)
sourceCpp('inst/causal/weighted_loglik.cpp')
## generate data stage 1
n <- 1000
p_1 <- 9
X_1 <- matrix(rnorm(n*p_1), nrow = n, ncol = p_1)
colnames(X_1) <- paste0('X1_', 1:p_1)

normal_log_prior <- function(theta){
  sum(dnorm(theta, 0, 1000, TRUE))
}

## generate binary treatment variable from logistic regression model given X_1
set.seed(1)
beta_true <- c(0.5, -1, rep(0, p_1-2))
lin_pred <- 1 + X_1 %*% beta_true
prob_treat <- 1/(1 + exp(-lin_pred))
hist(prob_treat)
treat <- rbinom(n, 1, prob_treat)

## prepare covariates for first regression
stage_1_model_mat <- cbind(1, X_1)
colnames(stage_1_model_mat)[1] <- 'intercept'
head(stage_1_model_mat)
p_1 <- ncol(stage_1_model_mat)
w1 <- rep(1, n)
## First posterior
logposterior1 <- function(theta1) normal_log_prior(theta1) + weighted_logistic_log_likC(theta1, treat, stage_1_model_mat, w1, p_1, n)
## Compute Laplace approximation
theta1hat <- (glm(treat ~ stage_1_model_mat -1, family="binomial")$coefficients) %>% unname
print(theta1hat)
theta1hat <- optim(par = theta1hat, fn = function(x) -logposterior1(x))$par
print(theta1hat)
## Compute Laplace approximation
asympvar1 <- solve(hessian(func = function(theta) -logposterior1(theta),
                           x = theta1hat))
laplace1_samples <- fast_rmvnorm(1e6, theta1hat, asympvar1)


## using the generic RWMRTH implementation
target <- list(dimension = p_1, density = function(theta1){
  return(sapply(1:dim(theta1)[1], function(ichain) logposterior1(theta1[ichain,])))
})

## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(rep(0.1, p_1))
tuning_parameters$niterations <- 5000
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
library(gridExtra)

df_stage1 <- (rbind(data.frame(method = "mcmc", theta = mcmc_samples_stage1),
                    data.frame(method = "laplace", theta = laplace1_samples[,1:p_1])))
g1 <- ggplot(df_stage1, aes(x = theta.1, fill = method)) +
  geom_density(alpha = 0.5)
g2 <- ggplot(df_stage1, aes(x = theta.2, fill = method)) +
  geom_density(alpha = 0.5)
g3 <- ggplot(df_stage1, aes(x = theta.3, fill = method)) +
  geom_density(alpha = 0.5)
grid.arrange(g1, g2, g3,  ncol = 3)

## generate outcome given treatment and propensity score
# introducing intercept
X_2 <- cbind(1, treat, lin_pred)
gamma_true <- c(1, 0.2, 0.1)
outcome <- X_2 %*% gamma_true + rnorm(n, 0, 1)
outcome <- outcome[,1]

N1 <- 1000
selected_indices <- round(seq(1, nrow(mcmc_samples_stage1), length.out = N1))
res_mcmc_list <- mclapply(selected_indices, function(i){
  beta_1 <- mcmc_samples_stage1[i,]
  prop_scores <- as.vector(stage_1_model_mat %*% matrix(beta_1, ncol = 1))
  mat_part_2 <- cbind(treat, prop_scores)
  p_2 <- ncol(mat_part_2) + 1
  ## calibrate proposal covariance from frequentist fit
  y = outcome
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
    rwmrth_res_2$chains[[ichain]][(tuning_parameters$niterations-9):tuning_parameters$niterations,]
  }
  return(mcmc_samples_stage2)
}, mc.cores = 10)

df_mcmc <- do.call('rbind', res_mcmc_list ) %>% as.data.frame()
df_mcmc$method <- 'cut Bayesian'
dim(df_mcmc)
head(df_mcmc)



logposterior2 <- function(theta1, theta2){
  prop_scores <- as.vector(stage_1_model_mat %*% matrix(theta1, ncol = 1))
  mat_part_2 <- cbind(treat, prop_scores)
  p_2 <- ncol(mat_part_2) + 1
  y = outcome
  x = mat_part_2
  prior_sd = 1000
  beta <- theta2[1:(p_2-1)]
  sigma <- theta2[p_2]
  if (sigma <= 0) return(-Inf)
  prior_sigma <- dgamma(sigma, 0.01, 0.01, log = TRUE)
  log_lik <- dnorm(y - as.vector(x%*%matrix(beta, ncol = 1)),
                   mean = 0, sd = sigma, log = TRUE) %>% sum()
  prior_beta <- sum(dnorm(beta, mean = 0, sd = prior_sd, log = TRUE))
  return(log_lik + prior_beta + prior_sigma)
}

theta1hat

theta2hat <- optim(par = c(rep(0,2),1),
                   fn = function(th2) -logposterior2(theta1hat, th2),
                   method = 'BFGS', control = list(reltol = 1e-10, maxit = 1000))$par

theta2hat

laplace_cut_res <- laplace_cut(logposterior1,
                                 logposterior2,
                                 init_theta1 = theta1hat,
                                 init_theta2 = theta2hat)
laplace_cut_samples <- fast_rmvnorm(1e6,
                                      laplace_cut_res$thetahat,
                                      laplace_cut_res$asympvar)


df_stage2 <- data.frame(method = "mcmc", theta = df_mcmc[,1:3] %>% as.matrix())
colnames(df_stage2) <- c("method", paste0("theta.", 1:3))
df_stage2 <- rbind(df_stage2, data.frame(method = "laplace", theta = laplace_cut_samples[,(p_1 + 1):(p_1 + 3)]))
g1 <- ggplot(df_stage2, aes(x = theta.1, fill = method)) +
  geom_density(alpha = 0.3) + xlab("theta 2 1")
g2 <- ggplot(df_stage2, aes(x = theta.2, fill = method)) +
  geom_density(alpha = 0.3) + xlab("theta 2 2")
g3 <- ggplot(df_stage2, aes(x = theta.3, fill = method)) +
  geom_density(alpha = 0.3) + xlab("theta 2 3")
grid.arrange(g1, g2, g3,  ncol = 3)
