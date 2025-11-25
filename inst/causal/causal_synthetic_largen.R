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

source('inst/causal/causal_inference_linear_regression_mcmc.R')
source('inst/causal/causal_inference_logistic_regression_mcmc.R')
## compile Rcpp function to evaluate (weighted) log-likelihood of logistic
## regression & linear regression
library(Rcpp)
sourceCpp('inst/causal/weighted_loglik.cpp')

## generate data stage 1
n <- 10000
p_1 <- 5
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


## generate outcome given treatment, propensity score and covariates
# introducing intercept
X_2 <- cbind(treat, lin_pred, X_1)
gamma_true <- c(0.2, 0.1, 0.5, -1, rep(0, p_1-3))
outcome <- X_2 %*% gamma_true + rnorm(n, 0, 1)
outcome <- outcome[,1]
# hist(outcome)

N1 <- 1000
selected_indices <- round(seq(1, nrow(mcmc_samples_stage1), length.out = N1))


nr_quantiles <- 5
res_mcmc_list <- mclapply(selected_indices, function(i){
  beta_1 <- mcmc_samples_stage1[i,]
  prop_scores_factor <- as.vector(stage_1_model_mat %*% matrix(beta_1, ncol = 1)) %>% gtools::quantcut(q = nr_quantiles)
  x_part_2 <- model.matrix(~prop_scores_factor) %>% as.matrix %>% unname()
  mat_part_2 <- cbind(treat, x_part_2)
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



### Two step Laplace approximation
# first stage
x <- stage_1_model_mat
p <- ncol(stage_1_model_mat)
n <- nrow(x)
theta1hat_glm <- (glm(treat ~ X_1, family="binomial")$coefficients) %>% unname
print(theta1hat_glm)
w1 <- rep(1, n)
logposterior1 <- function(theta1) normal_log_prior(theta1) + weighted_logistic_log_likC(theta1, treat, stage_1_model_mat, w1, p_1, n)
theta1hat <- optim(par = theta1hat_glm,
                   fn = function(x) -logposterior1(x))$par
print(theta1hat)

# second stage
prop_scores_factor <- as.vector(stage_1_model_mat%*%theta1hat) %>% gtools::quantcut(q = nr_quantiles)
x_part_2 <- model.matrix(~prop_scores_factor) %>%
  as.data.frame() %>%
  as.matrix %>%
  unname()

mat_part_2 <- cbind(treat, x_part_2)
p_2 <- ncol(mat_part_2) + 1
## calibrate proposal covariance from frequentist fit
y = outcome
x = mat_part_2
prior_sd = 1000
freq_model <- lm(y~x-1)
betahat <- freq_model %>% coefficients
sigmahat <- sqrt(freq_model$residuals %>% var())

log_target <- function(theta2){
  beta <- theta2[1:(p_2-1)]
  sigma <- theta2[p_2]
  if (sigma <= 0) return(-Inf)
  prior_sigma <- dgamma(sigma, 0.01, 0.01, log = TRUE)
  log_lik <- dnorm(y - as.vector(x%*%matrix(beta, ncol = 1)),
                   mean = 0, sd = sigma, log = TRUE) %>% sum()
  prior_beta <- sum(dnorm(beta, mean = 0, sd = prior_sd, log = TRUE))
  return(log_lik + prior_beta + prior_sigma)
}

start_theta2 <- c(betahat, sigmahat)
theta2hat <- optim(par = start_theta2,
      fn = function(th2) -log_target(th2),
      method = 'BFGS', control = list(reltol = 1e-10, maxit = 1000))$par

theta2hat


logposterior2 <- function(theta1, theta2){
  prop_scores_factor <- as.vector(stage_1_model_mat %*% matrix(theta1, ncol = 1)) %>% gtools::quantcut(q = nr_quantiles)
  x_part_2 <- model.matrix(~prop_scores_factor) %>%
    as.matrix %>%
    unname()
  mat_part_2 <- cbind(treat, x_part_2)
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

d1 <- p_1
d2 <- p_2
xseq <- seq(from = -20, to = -10, by = 0.1)
yseq <- sapply(xseq,
       function(x) logposterior2(theta1hat+ exp(x) * rep(1, d1), theta2hat))
plot(xseq, yseq, type = "l")

## Computes hessians at the two stage MAP estimate
hessian1 <- hessian(func = logposterior1, theta1hat)
## Compute Hessian of L_2 at 2SMAP
hessian2 <- hessian(func = function(theta) logposterior2(theta[1:d1], theta[(d1+1):(d1+d2)]), c(theta1hat, theta2hat),
                    method.args=list(eps=1e-1, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=8, v=2, show.details=FALSE))

numDeriv::grad(func = function(th) logposterior2(th[1:d1], th[(d1+1):(d1+d2)]), x = c(theta1hat,theta2hat),
               method.args=list(eps=1e-2, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=10, v=2, show.details=FALSE))


## Construct asymptotic precision matrix of the Laplace approximation
asympprec <- matrix(0, nrow = d1 + d2, ncol = d1 + d2)
topleft <- hessian1 + hessian2[1:d1,(d1+1):(d1+d2),drop=F] %*% solve(hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2),drop=F], hessian2[(d1+1):(d1+d2),1:d1,drop=F])
asympprec[1:d1, 1:d1] <- topleft
asympprec[1:d1,(d1+1):(d1+d2)] <- hessian2[1:d1,(d1+1):(d1+d2)]
asympprec[(d1+1):(d1+d2),1:d1] <- hessian2[(d1+1):(d1+d2),1:d1]
asympprec[(d1+1):(d1+d2),(d1+1):(d1+d2)] <- hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2)]
## Laplace approximation of cut posterior is given by N(thetahat, asympvar)
laplace_cut_asympvar <- solve(-asympprec)
laplace_cut_asympvar[(p_1+1):(p_1+2),(p_1+1):(p_1+2)]
cov(df_mcmc[,1:2])

## alternate
asympvar_2 <- matrix(0, nrow = d1 + d2, ncol = d1 + d2)
asympvar_2[1:d1, 1:d1] <- solve(hessian1)
# laplace_cut_asympvar[1:d1, 1:d1]
asympvar_2[1:d1, (d1+1):(d1+d2)] <-  - solve(hessian1) %*% hessian2[1:d1,(d1+1):(d1+d2),drop=F] %*% solve(hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2),drop=F])
# asympvar_2[1:d1, (d1+1):(d1+d2)]
# laplace_cut_asympvar[1:d1, (d1+1):(d1+d2)]
asympvar_2[(d1+1):(d1+d2),1:d1] <-  t(asympvar_2[1:d1, (d1+1):(d1+d2)])
# asympvar_2[(d1+1):(d1+d2),1:d1]
# laplace_cut_asympvar[(d1+1):(d1+d2),1:d1]
asympprec[(d1+1):(d1+d2),(d1+1):(d1+d2)] <- solve(hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2),drop=F]) - asympvar_2[(d1+1):(d1+d2),1:d1] %*% hessian2[1:d1,(d1+1):(d1+d2),drop=F] %*% solve(hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2),drop=F])
laplace_cut_asympvar[(d1+1):(d1+d2),(d1+1):(d1+d2)]

laplace_cut_results <- laplace_cut(logposterior1, logposterior2, theta1hat, theta2hat)
laplace_cut_results$asympvar[(p_1+1):(p_1+2),(p_1+1):(p_1+2)]
cov(df_mcmc[,1:7])

as.numeric(laplace_cut_results$thetahat)
as.numeric(c(colMeans(mcmc_samples_stage1), colMeans(df_mcmc[,1:p_2])))

save(laplace_cut_results, theta1hat, theta2hat, mcmc_samples_stage1, df_mcmc,
     file = 'inst/causal/causal_synthetic_largen.RData')

load(file = 'inst/causal/causal_synthetic_largen.RData')

## generate samples from cut-Laplace
laplace_cut_samples <- fast_rmvnorm(1e5, c(theta1hat, theta2hat), laplace_cut_results$asympvar)
head(laplace_cut_samples)
head(df_mcmc)

## first stage accuracy
library(gridExtra)

df_stage1 <- (rbind(data.frame(method = "mcmc", theta = mcmc_samples_stage1),
      data.frame(method = "laplace", theta = laplace_cut_samples[,1:p_1]),
      data.frame(method = "laplace1", theta = laplace1_samples[,1:p_1])))
g1 <- ggplot(df_stage1, aes(x = theta.1, fill = method)) +
  geom_density(alpha = 0.3)
g2 <- ggplot(df_stage1, aes(x = theta.2, fill = method)) +
  geom_density(alpha = 0.3)
g3 <- ggplot(df_stage1, aes(x = theta.3, fill = method)) +
  geom_density(alpha = 0.3)
grid.arrange(g1, g2, g3,  ncol = 3)

## second stage accuracy
## Laplace approximation at second stage, for fixed theta1hat
Laplace2giventheta1hat <- solve(-hessian(func = function(th2) logposterior2(theta1hat, th2),
                                         x = theta2hat))
Laplace2giventheta1hat[1:2,1:2]
cov(df_mcmc[,1:2])


df_stage2 <- data.frame(method = "mcmc", theta = df_mcmc[,1:6] %>% as.matrix())
colnames(df_stage2) <- c("method", paste0("theta.", 1:6))
df_stage2 <- rbind(df_stage2, data.frame(method = "laplace", theta = laplace_cut_samples[,(p_1 + 1):(p_1 + 6)]))
laplace2_samples <- fast_rmvnorm(1e5, theta2hat, Laplace2giventheta1hat)
df_stage2 <- rbind(df_stage2, data.frame(method = "laplace2", theta = laplace2_samples[,1:6]))
g1 <- ggplot(df_stage2, aes(x = theta.1, fill = method)) +
  geom_density(alpha = 0.3) + xlab("theta 2 1")
g2 <- ggplot(df_stage2, aes(x = theta.2, fill = method)) +
  geom_density(alpha = 0.3) + xlab("theta 2 2")
g3 <- ggplot(df_stage2, aes(x = theta.3, fill = method)) +
  geom_density(alpha = 0.3) + xlab("theta 2 3")
grid.arrange(g1, g2, g3,  ncol = 3)
# g4 <- ggplot(df_stage2, aes(x = theta.4, fill = method)) +
#   geom_density(alpha = 0.7)
# g5 <- ggplot(df_stage2, aes(x = theta.5, fill = method)) +
#   geom_density(alpha = 0.7)
# g6 <- ggplot(df_stage2, aes(x = theta.6, fill = method)) +
#   geom_density(alpha = 0.7)
# grid.arrange(g1, g2, g3, g4, g5, g6,  ncol = 3)


df_compar <- data.frame(V1 = laplace_cut_samples[,p_1 + 1], method = 'cut Laplace')
df_compar <- rbind(df_compar, df_mcmc %>% dplyr::select(V1, method))
## histogram of 'treatment' coefficient
gtreat_largen <- ggplot(df_compar, aes(x = V1, fill = method)) +
  geom_density(alpha = 0.7) +
  labs(x = 'Treatment effect',
       y = 'Density')
gtreat_largen
ggsave('inst/causal/causal_treat_largen.pdf', gtreat_largen, width = 8, height = 5)

