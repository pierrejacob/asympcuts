## load packages
library(asympcuts)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(MatchIt)
library(BayesLogit) # Polya-Gamma sampler rpg()
## load LaplacesDemon for the inverse-Gamma distribution
library(LaplacesDemon)

## parallel-safe random number streams for mclapply
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

## compile Rcpp function to evaluate (weighted) log-likelihood of logistic
## regression & linear regression
## - weighted_logistic_log_likC(beta, y, x, weights, ncol, nrow)
## - weighted_linear_regressionC(beta, y, x, weights, ncol, nrow)
library(Rcpp)
sourceCpp('inst/causal/causal_lalonde/weighted_loglik.cpp')

#### Data preparation ####
## load data 'lalonde' from MatchIt
## Description: this is a subsample of the data from the treated group in the
## National Supported Work Demonstration (NSW) and the comparison sample from
## the Population Survey of Income Dynamics (PSID). This data was previously
## analyzed extensively by Lalonde (1986) and Dehejia and Wahba (1999).
data(lalonde)

## Standardize the data
## Can we use the 'scale' function in the base package?
standardize <- function(x) {
  (x-mean(x))/sd(x)
}

## have a look
head(lalonde)
# preparation
lalonde$race <- factor(lalonde$race, levels = c('white', 'hispan', 'black'))
lalonde$age <- standardize(lalonde$age)
lalonde$educ <- standardize(lalonde$educ)
lalonde$re74 <- standardize(lalonde$re74)
lalonde$re75 <- standardize(lalonde$re75)
lalonde$age_sq <- standardize((lalonde$age)^2)

#### Stage 1: logistic regression - MAP estimation ####

# Normal(0,100^2) prior for all coefficients in model 1
prior1_sd <- 100
log_prior1 <- function(theta){
  sum(dnorm(theta, 0, prior1_sd, TRUE))
}
## Outcome: `treat`
## Covariates: age, age_sq, educ, race, married, nodegree, re74, re75 (8 covariates)
## the model.matrix function adds an intercept and converts race into two indicators (racehispan, raceblack)
stage_1_model_mat <- model.matrix(treat ~ age + age_sq + educ + race + married + nodegree + re74 + re75,
                                  data = lalonde)
head(stage_1_model_mat)
p_1 <- ncol(stage_1_model_mat)
cat("Number of covariates in stage 1: ", p_1, "\n")
n <- nrow(stage_1_model_mat)
cat("Number of observations in stage 1: ", n, "\n")
w1 <- rep(1, n)
# compute log-posterior density
logposterior1 <- function(theta1) log_prior1(theta1) + weighted_logistic_log_likC(theta1, lalonde$treat, stage_1_model_mat, w1, p_1, n)
## compute MLE
theta1mle <- (glm(lalonde$treat ~ stage_1_model_mat -1, family="binomial")$coefficients) %>% unname
## compute MAP estimate, only slightly different
theta1hat <- optim(par = theta1mle, fn = function(x) -logposterior1(x), control = list(reltol = 1e-10))$par
print(theta1hat) # MAP estimate

## Compute variance of the Laplace approximation
asympvar1 <- solve(hessian(func = function(theta) -logposterior1(theta),
                           x = theta1hat))
## The Laplace approximation is N(theta1hat, asympvar1)

#### Cut Bayesian approach ####

#### Stage 1: logistic regression - Polya-Gamma Gibbs sampler ####
## Polson, Scott & Windle (2013, JASA): with kappa_i = y_i - 1/2 and
## latent omega_i ~ PG(1, x_i' beta), the full conditionals are
##   omega_i | beta    ~ PG(1, x_i' beta)
##   beta    | omega   ~ N(m, V),  V = (X' Omega X + B^{-1})^{-1},  m = V X' kappa
## with prior beta ~ N(0, B), here B = prior1_sd^2 I (same prior as above).
## Polya-Gamma draws from BayesLogit::rpg (install.packages("BayesLogit")).
pg_gibbs_logistic <- function(y, X, niterations, beta_init, prior_sd = prior1_sd){
  p <- ncol(X); n <- nrow(X)
  Xtkappa <- crossprod(X, y - 0.5)
  prior_prec <- diag(1 / prior_sd^2, p)
  beta <- beta_init
  chain <- matrix(NA_real_, nrow = niterations, ncol = p)
  chain[1, ] <- beta
  for (it in 2:niterations){
    ## omega | beta
    omega <- BayesLogit::rpg(num = n, h = 1, z = as.vector(X %*% beta))
    ## beta | omega: posterior precision Q = X' Omega X + B^{-1} = R'R
    R <- chol(crossprod(X * sqrt(omega)) + prior_prec)
    mu <- backsolve(R, forwardsolve(t(R), Xtkappa))
    beta <- as.vector(mu + backsolve(R, rnorm(p)))
    chain[it, ] <- beta
  }
  return(chain)
}

## tuning parameters
tuning_parameters <- list()
tuning_parameters$niterations <- 10000
tuning_parameters$nchains <- 4
tuning_parameters$rinit <- function(nchains){
  chains <- matrix(rnorm(nchains*p_1, 0, 1), ncol = p_1)
  return(chains)
}
## initialization
init_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
## Gibbs runs, one per chain
pg_res_1 <- lapply(1:tuning_parameters$nchains, function(ichain)
  pg_gibbs_logistic(lalonde$treat, stage_1_model_mat,
                    niterations = tuning_parameters$niterations,
                    beta_init = init_chains[ichain, ]))
## conservative choice of burnin (PG Gibbs typically converges within a few iterations)
burnin <- tuning_parameters$niterations/2

mcmc_samples_stage1 <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
  pg_res_1[[ichain]][(burnin+1):tuning_parameters$niterations,]
}
## quick checks: posterior mean vs MAP, and trace plot of the first chain
# print(rbind(pg_mean = colMeans(mcmc_samples_stage1), map = theta1hat))
# matplot(pg_res_1[[1]], type = 'l')
## MCMC samples approximating first posterior
print(dim(mcmc_samples_stage1))

## first stage: compare MCMC samples and Laplace approximation
laplace1_samples <- fast_rmvnorm(1e6, theta1hat, asympvar1)

df_stage1 <- (rbind(data.frame(method = "mcmc", theta = mcmc_samples_stage1),
                    data.frame(method = "laplace", theta = laplace1_samples[,1:p_1])))

## reshape and plot marginals of each theta.X overlaid by method
df_stage1_forplotting <- df_stage1 %>%
  tidyr::pivot_longer(cols = starts_with("theta"), names_to = "theta", values_to = "value") %>%
  dplyr::mutate(theta = factor(theta, levels = paste0("theta.", 1:p_1)))

# g_mcmc_stage1 <- ggplot(df_stage1_forplotting, aes(x = value, fill = method)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~theta, scales = "free") +
#   labs(x = "Value", y = "Density", fill = "Method") +
#   theme_minimal()
# print(g_mcmc_stage1)
## as expected with a large sample size,
## the Laplace approximation is very close to the MCMC samples

## select N1 samples from the MCMC samples of stage 1 to run stage 2
N1 <- 5000
selected_indices <- round(seq(1, nrow(mcmc_samples_stage1), length.out = N1))

#### Stage 2: linear regression - 2SMAP ####
sd_re78 <- sd(lalonde$re78)
lalonde$re78 <- standardize(lalonde$re78)
## Covariates?
## Strata built in the spirit of Dehejia & Wahba (1999, see ~page 8):
## 1. common support: discard the controls whose estimated score is below the
##    smallest score among the treated (no treated unit is ever discarded);
## 2. stratum boundaries = quantiles of the score among the treated units, so
##    that each stratum contains about n_T / nqtles treated units;
## 3. a stratum without controls (or without treated) is merged with a
##    neighbouring stratum, so the number of strata K can be < nqtles.
## The linear predictor is used instead of the score itself: it is a monotone
## transformation of the score, so the trimming and the strata are the same.
## (D&W also refine the strata until covariates are balanced within each
## stratum; that step is not reproduced here for simplicity.)
## Prior: Normal(0,100^2) on each parameter + Inverse-Gamma(0.1,0.1) on the sigma2
prior2_coef_sd <- 100 # on each regression coefficient
prior2_invgamma_param1  <- 1
prior2_invgamma_param2  <- 1
treatment <- lalonde$treat
## create matrix of covariates from theta1; needs treatment variable in global environment
## Returns the model matrix (K stratum indicators + K interactions with the
## treatment, no intercept) for the retained units only; the row indices of the
## retained units in lalonde are stored in attr(, "keep").
create_covariate2 <- function(theta1, nqtles = 5){
  lin_pred <- as.vector(stage_1_model_mat %*% theta1)
  ## 1. common support
  keep <- which(lin_pred >= min(lin_pred[treatment == 1]) | lin_pred <= max(lin_pred[treatment == 1]))

  lp <- lin_pred[keep]
  tr <- treatment[keep]
  ## 2. interior boundaries = quantiles of the score among the treated;
  ## strata are (b_{k-1}, b_k], the lowest one also contains the minimum
  breaks <- quantile(lp[tr == 1], probs = seq(0, 1, length.out = nqtles + 1), names = FALSE)
  breaks <- unique(breaks[-c(1, nqtles + 1)])
  stratum <- findInterval(lp, breaks, left.open = TRUE) + 1
  ## 3. merge strata without controls or without treated with a neighbour
  repeat {
    stratum <- match(stratum, sort(unique(stratum))) # relabel as 1..K
    K <- max(stratum)
    ncontrols <- tabulate(stratum[tr == 0], nbins = K)
    ntreated <- tabulate(stratum[tr == 1], nbins = K)
    bad <- which(ncontrols == 0 | ntreated == 0)
    if (length(bad) == 0 || K == 1) break
    k <- bad[1]
    stratum[stratum == k] <- if (k > 1) k - 1 else 2
  }
  ## model matrix: stratum indicators and their interactions with treatment
  S <- outer(stratum, 1:K, "==") * 1
  covariates_2 <- cbind(S, S * tr)
  colnames(covariates_2) <- c(paste0("S", 1:K), paste0("S", 1:K, "xT"))
  attr(covariates_2, "keep") <- keep
  return(covariates_2)
}

## ATT = sum_k (n_k^T / n_T) * (treated - control difference in stratum k)
## theta2: vector whose entries K+1, ..., 2K are the interaction coefficients
## (extra entries, e.g. sigma, are ignored)
att_from_beta <- function(theta2, covariates_2){
  K <- ncol(covariates_2) / 2
  idx <- K + seq_len(K)
  ntreated_k <- colSums(covariates_2[, idx, drop = FALSE])
  return(sum(ntreated_k / sum(ntreated_k) * theta2[idx]))
}

covariates_fromtheta1hat <- create_covariate2(theta1hat)
keep_hat <- attr(covariates_fromtheta1hat, "keep")
y_hat <- lalonde$re78[keep_hat]
K_hat <- ncol(covariates_fromtheta1hat) / 2
cat("Controls kept at theta1hat:", sum(treatment[keep_hat] == 0), "out of", sum(treatment == 0), "\n")
## number of treated and controls per stratum at theta1hat
ntreated_hat <- colSums(covariates_fromtheta1hat[, K_hat + 1:K_hat, drop = FALSE])
print(rbind(treated = ntreated_hat,
            controls = colSums(covariates_fromtheta1hat[, 1:K_hat, drop = FALSE]) - ntreated_hat))

# log posterior for the second stage, given covariates
# argument: vector of p_2 regression coefficients + std deviation of residuals
# y defaults to the outcomes of the units retained by create_covariate2

lp2_ <- function(theta2, covariates_2, y = lalonde$re78[attr(covariates_2, "keep")]){
  p_2 <- ncol(covariates_2)
  n_2 <- nrow(covariates_2)
  beta2 <- theta2[1:p_2]
  sigma_ <- theta2[p_2+1]
  ## prior is Normal(0,prior2_coef_sd^2) on beta[j]
  ## and Inverse-Gamma(prior2_invgamma_param1,prior2_invgamma_param2) on sigma2
  if (sigma_ <= 0) return(-Inf)
  ## note that the prior is placed on the variance of the residuals
  prior_part <- sum(dnorm(beta2, mean = 0, sd = prior2_coef_sd, log = TRUE)) +
    dinvgamma(sigma_^2, prior2_invgamma_param1, prior2_invgamma_param2, log = TRUE)
  ll_part <- - n_2 * log(sigma_) +
    (1/sigma_^2) * weighted_linear_regressionC(theta2, y, covariates_2, rep(1, n_2), p_2, n_2)
  return(prior_part + ll_part)
}
## two-stage MAP
p_2hat <- ncol(covariates_fromtheta1hat)
theta2hat <- optim(par = c(rep(0, p_2hat), 1), fn = function(x) -lp2_(x, covariates_fromtheta1hat),
                   method = "BFGS", control = list(reltol = 1e-12, maxit = 1000))$par
## compare with OLS
ols_hat <- lm(y_hat ~ covariates_fromtheta1hat - 1)
print(rbind(map = theta2hat[1:p_2hat], ols = coef(ols_hat)))
## ATT using 2SMAP
att_hat <- att_from_beta(theta2hat, covariates_fromtheta1hat)
cat("ATT using two-stage MAP:", att_hat, "standardized,", att_hat * sd_re78, "dollars\n")

#### Stage 2: exact Gibbs sampler for the linear regression ####
## Model:  y ~ N(x beta, sigma2 I)
## Prior:  beta ~ N(0, prior2_sd^2 I),  sigma2 ~ InvGamma(a0, b0)
## NB: this differs from lp2_ (which puts a Gamma(0.01, 0.01) prior on sigma);
## the inverse-gamma prior on sigma2 makes both full conditionals available
## in closed form. With n = 614 the choice of prior is immaterial.
## Columns of zeros (strata without treated units) need no special treatment:
## the prior precision keeps the posterior precision matrix invertible, and the
## corresponding coefficients are drawn from their prior.
gibbs_stage2 <- function(y, x, niter = 50, nkeep = 2,
                         prior_sd = prior2_coef_sd, a0 = 1, b0 = 1){
  p <- ncol(x); n <- nrow(x)
  XtX <- crossprod(x)
  Xty <- crossprod(x, y)
  prior_prec <- diag(1 / prior_sd^2, p)
  sigma2 <- var(y) # initial value
  draws <- matrix(NA_real_, nrow = nkeep, ncol = p + 1)
  for (it in 1:niter){
    ## beta | sigma2, y ~ N(mu, Q^{-1}) with Q = X'X/sigma2 + prior precision
    R <- chol(XtX / sigma2 + prior_prec) # Q = R'R
    mu <- backsolve(R, forwardsolve(t(R), Xty / sigma2))
    beta <- as.vector(mu + backsolve(R, rnorm(p)))
    ## sigma2 | beta, y ~ InvGamma(a0 + n/2, b0 + RSS/2)
    rss <- sum((y - as.vector(x %*% beta))^2)
    sigma2 <- 1 / rgamma(1, shape = a0 + n / 2, rate = b0 + rss / 2)
    if (it > niter - nkeep) draws[it - (niter - nkeep), ] <- c(beta, sqrt(sigma2))
  }
  return(draws)
}

## sanity check at theta1hat: Gibbs posterior mean vs OLS
# gibbs_check <- gibbs_stage2(y_hat, covariates_fromtheta1hat, niter = 5000, nkeep = 4000)
# rbind(gibbs = colMeans(gibbs_check), ols = c(coef(ols_hat), NA))

## MCMC for stage 2 using fixed theta1hat
# covariates_fromtheta1hat <- create_covariate2(theta1hat)
# keep_hat <- attr(covariates_fromtheta1hat, "keep")
# y_hat <- lalonde$re78[keep_hat]
# gibbs_stage2_theta1hat <- gibbs_stage2(y_hat, covariates_fromtheta1hat,
#              niter = 1e4, nkeep = 5e3, prior_sd = prior2_coef_sd,
#              a0 = prior2_invgamma_param1, b0 = prior2_invgamma_param2)
# head(gibbs_stage2_theta1hat)
# att_theta1hat <- apply(gibbs_stage2_theta1hat, 1, function(theta2) att_from_beta(theta2, covariates_fromtheta1hat))
# head(att_theta1hat)

## number of iterations per theta1 and number of draws kept (the last ones)
niter_gibbs <- 50
nkeep_gibbs <- 2

res_mcmc_list <- mclapply(selected_indices, function(i){
  beta_1 <- mcmc_samples_stage1[i,]
  ## create the matrix of covariates for the second stage (retained units only)
  covariates_2 <- create_covariate2(beta_1)
  keep <- attr(covariates_2, "keep")
  y_2 <- lalonde$re78[keep]
  ## Gibbs sampler for the second-stage posterior given theta1
  mcmc_samples_stage2 <- gibbs_stage2(y_2, covariates_2,
                                    niter = niter_gibbs, nkeep = nkeep_gibbs, prior_sd = prior2_coef_sd,
                                    a0 = prior2_invgamma_param1, b0 = prior2_invgamma_param2)
  ## Average treatment on the treated (ATT):
  ## sum_{k=1}^K (n_k^T/n_T) * (beta_{K+k})
  att <- apply(mcmc_samples_stage2, 1, function(theta2) att_from_beta(theta2, covariates_2))
  return(list(att = matrix(att, ncol = 1),
              diagnostics = c(nstrata = ncol(covariates_2) / 2,
                              ncontrols_kept = sum(treatment[keep] == 0))))
}, mc.cores = 10)

## diagnostics: number of strata after merging, number of controls retained
df_diag <- do.call('rbind', lapply(res_mcmc_list, function(r) r$diagnostics)) %>% as.data.frame()
print(table(df_diag$nstrata))
print(summary(df_diag$ncontrols_kept))

df_mcmc <- do.call('rbind', lapply(res_mcmc_list, function(r) r$att)) %>% as.data.frame()
df_mcmc$method <- "cut Bayesian (PG Gibbs + Gibbs, D&W strata)"

tail(df_mcmc)
hist(df_mcmc$V1, nclass = 100)

## average treatment effect on the treated (ATT) is in the original scale of re78
hist(df_mcmc$V1 * sd_re78, nclass = 100)
##
mean(df_mcmc$V1 * sd_re78)


#### Posterior Bootstrap for Modular Inference (PBMI) ####
## Each PBMI draw is obtained as follows.
## Stage 1: draw weights w_i ~ Exp(1), i = 1, ..., n, and compute
##   theta1 = argmax  sum_i w_i log p(treat_i | x_i, theta1) + w_0 * log prior1(theta1)
## Stage 2: build the strata from theta1 with create_covariate2, which retains
##   the units in `keep` (common support), and compute
##   theta2 = argmax  sum_{i in keep} v_i log N(re78_i; x2_i' beta2, sigma^2)
##                    + v_0 * log prior2(theta2)
##   where the stage 2 weights are
##   - "same weights":       v_i = w_i for i in keep, i.e. the retained
##                            individuals keep their stage 1 weights;
##   - "refreshed weights":  v_i ~ Exp(1), drawn independently of the w_i.
## Then the ATT is computed from theta2 exactly as in the cut Bayesian approach
## (att_from_beta, with stratum sizes n_k^T computed without weights).
##
## w_0 and v_0 scale the log-priors. With w_0 = v_0 = 1 the targets are the
## log-posteriors used for 2SMAP above, so that with all weights equal to 1 the
## procedure returns (theta1hat, theta2hat). Setting w_0 = v_0 = 0 ignores the priors.
pb_w0 <- 1
pb_v0 <- 1
## number of PBMI draws: as many as the number of cut Bayesian draws
N_pb <- nrow(df_mcmc)
pb_nr_cores <- 10
pb_control <- list(reltol = 1e-12, maxit = 1000)

## Stage 1: negative weighted log-target and its gradient
pb_negtarget1 <- function(theta1, w, w_0 = pb_w0){
  - weighted_logistic_log_likC(theta1, lalonde$treat, stage_1_model_mat, w, p_1, n) -
    w_0 * log_prior1(theta1)
}
pb_neggrad1 <- function(theta1, w, w_0 = pb_w0){
  prob <- plogis(as.vector(stage_1_model_mat %*% theta1))
  - as.vector(crossprod(stage_1_model_mat, w * (lalonde$treat - prob))) +
    w_0 * theta1 / prior1_sd^2
}

## Stage 2: negative weighted log-target and its gradient, parametrized by
## par = (beta2, log(sigma)) so that the optimization is unconstrained.
## The function of sigma is the same as in lp2_ (no Jacobian term), so the
## maximizer in sigma is the same as with the sigma parametrization.
## With v = 1 and v_0 = 1 this equals -lp2_(c(beta2, sigma), covariates_2, y).
pb_negtarget2 <- function(par, covariates_2, y, v, v_0 = pb_v0){
  p_2 <- ncol(covariates_2)
  beta2 <- par[1:p_2]
  sigma_ <- exp(par[p_2 + 1])
  resid <- y - as.vector(covariates_2 %*% beta2)
  loglik <- - sum(v) * log(sigma_) - 0.5 * sum(v * resid^2) / sigma_^2
  logprior <- sum(dnorm(beta2, mean = 0, sd = prior2_coef_sd, log = TRUE)) +
    dinvgamma(sigma_^2, prior2_invgamma_param1, prior2_invgamma_param2, log = TRUE)
  return(- loglik - v_0 * logprior)
}
pb_neggrad2 <- function(par, covariates_2, y, v, v_0 = pb_v0){
  p_2 <- ncol(covariates_2)
  beta2 <- par[1:p_2]
  sigma_ <- exp(par[p_2 + 1])
  resid <- y - as.vector(covariates_2 %*% beta2)
  grad_beta <- as.vector(crossprod(covariates_2, v * resid)) / sigma_^2 -
    v_0 * beta2 / prior2_coef_sd^2
  ## derivative with respect to log(sigma); for the prior,
  ## log dinvgamma(sigma^2; a, b) = const - (a + 1) * 2 * log(sigma) - b / sigma^2
  grad_logsigma <- - sum(v) + sum(v * resid^2) / sigma_^2 +
    v_0 * (- 2 * (prior2_invgamma_param1 + 1) + 2 * prior2_invgamma_param2 / sigma_^2)
  return(- c(grad_beta, grad_logsigma))
}

## One PBMI draw. Returns theta1, theta2 (= beta2, sigma), the ATT and diagnostics;
## theta2 and the ATT are NA if one of the optimization algorithms did not converge.
pbmi_draw <- function(same_weights, w = rexp(n, 1)){
  ## stage 1, started from the MAP
  opt1 <- optim(par = theta1hat, fn = pb_negtarget1, gr = pb_neggrad1, w = w,
                method = "BFGS", control = pb_control)
  theta1 <- opt1$par
  ## strata and retained units given theta1
  covariates_2 <- create_covariate2(theta1)
  keep <- attr(covariates_2, "keep")
  y_2 <- lalonde$re78[keep]
  p_2 <- ncol(covariates_2)
  ## stage 2 weights
  if (same_weights){
    v <- w[keep]
  } else {
    v <- rexp(length(keep), 1)
  }
  ## stage 2, started from weighted least squares
  wls <- lm.wfit(covariates_2, y_2, v)
  par_init <- c(wls$coefficients, log(sqrt(sum(v * wls$residuals^2) / sum(v))))
  opt2 <- optim(par = par_init, fn = pb_negtarget2, gr = pb_neggrad2,
                covariates_2 = covariates_2, y = y_2, v = v,
                method = "BFGS", control = pb_control)
  theta2 <- c(opt2$par[1:p_2], exp(opt2$par[p_2 + 1]))
  converged <- (opt1$convergence == 0) && (opt2$convergence == 0)
  att <- if (converged) att_from_beta(theta2, covariates_2) else NA_real_
  return(list(theta1 = theta1,
              theta2 = theta2,
              att = att,
              diagnostics = c(nstrata = p_2 / 2,
                              ncontrols_kept = sum(treatment[keep] == 0),
                              converged = converged)))
}

## sanity check: with all weights equal to one (and w_0 = v_0 = 1),
## PBMI returns the two-stage MAP (theta1hat, theta2hat)
if (pb_w0 == 1 && pb_v0 == 1){
  check_unit_weights <- pbmi_draw(same_weights = TRUE, w = rep(1, n))
  cat("Unit weights vs 2SMAP, max abs difference: theta1",
      max(abs(check_unit_weights$theta1 - theta1hat)),
      ", theta2", max(abs(check_unit_weights$theta2 - theta2hat)), "\n")
}

## run PBMI; returns the stage 1 draws, the ATT draws and diagnostics
run_pbmi <- function(N, same_weights, nr_cores = pb_nr_cores){
  res <- mclapply(1:N, function(iter) pbmi_draw(same_weights), mc.cores = nr_cores)
  list(theta1 = do.call('rbind', lapply(res, function(r) r$theta1)),
       att = sapply(res, function(r) r$att),
       diagnostics = do.call('rbind', lapply(res, function(r) r$diagnostics)) %>% as.data.frame())
}

## same weights: retained individuals keep their stage 1 weights
pb_res <- run_pbmi(N_pb, same_weights = TRUE)
## independently refreshed weights in stage 2
pb_res_refresh <- run_pbmi(N_pb, same_weights = FALSE)

## diagnostics
for (pbres in list(pb_res, pb_res_refresh)){
  cat("Non-converged draws:", sum(pbres$diagnostics$converged == 0), "out of", nrow(pbres$diagnostics), "\n")
  print(table(pbres$diagnostics$nstrata))
  print(summary(pbres$diagnostics$ncontrols_kept))
}

## stage 1: compare PBMI draws with the PG Gibbs draws (means and sds)
print(rbind(pg_mean = colMeans(mcmc_samples_stage1), pb_mean = colMeans(pb_res$theta1),
            pg_sd = apply(mcmc_samples_stage1, 2, sd), pb_sd = apply(pb_res$theta1, 2, sd)))

df_pb <- data.frame(V1 = pb_res$att, method = "PBMI (same weights)")
df_pb_refresh <- data.frame(V1 = pb_res_refresh$att, method = "PBMI (refreshed weights)")

# df_stage2giventheta1hat <- data.frame(V1 = att_theta1hat, method = "Stage2giventheta1hat")

df_twostage <- rbind(df_mcmc, df_pb, df_pb_refresh) %>% dplyr::filter(!is.na(V1))

## ATT in dollars, by method
df_twostage %>% group_by(method) %>%
  summarise(mean = mean(V1 * sd_re78), sd = sd(V1 * sd_re78),
            q2.5 = quantile(V1 * sd_re78, 0.025), q97.5 = quantile(V1 * sd_re78, 0.975),
            n = n()) %>%
  print()


## Visualization to compare PBMI and cut Bayesian
g_att <- ggplot(df_twostage, aes(x = V1 * sd_re78, colour = method)) +
  geom_density() +
  geom_vline(xintercept = att_hat * sd_re78, linetype = 2) +
  labs(x = "ATT (dollars)", y = "Density", colour = "Method") +
  theme_minimal()
print(g_att)

ggsave('inst/causal/causal_lalonde/causal_att.pdf', g_att, width = 8, height = 5)

save(theta1hat, theta2hat, att_hat, sd_re78, df_twostage, pb_res, pb_res_refresh,
     file = 'inst/causal/causal_lalonde/causal_twostage_v3_pbmi_results.RData')
