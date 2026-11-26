# Logistic propensity posterior, sampled by Stan through rstanarm.
# ndpost is the TOTAL retained draws; nskip is warmup PER CHAIN.
fit_logistic_posterior <- function(z, a, ndpost, nskip, chains, seed) {
  if (!requireNamespace("rstanarm", quietly = TRUE))
    stop("Install rstanarm first: install.packages('rstanarm')")
  controls <- c(ndpost, nskip, chains, seed)
  stopifnot(all(is.finite(controls)), all(controls == floor(controls)),
            chains >= 1, nskip >= 0, ndpost %% chains == 0,
            ndpost / chains >= 20, seed > 0, seed <= .Machine$integer.max)
  stopifnot(is.matrix(z), is.numeric(z), ncol(z) >= 1,
            all(is.finite(z)), all(z[, 1L] == 1),
            length(a) == nrow(z), all(a %in% c(0, 1)))
  dimension <- ncol(z)
  retained <- ndpost %/% chains # calculate retained draws per chain
  colnames(z) <- paste0("alpha_", 0:(dimension - 1L))

  # The first column is already the intercept. `0 +` prevents another one
  # and preserves our prior on the intercept WITHOUT centering predictors.
  # Priors are independent N(0, 10^2) for alpha_0 and N(0, 5^2) otherwise.
  fit <- rstanarm::stan_glm(
    a ~ 0 + ., data = data.frame(a = as.integer(a), z),
    family = binomial(link = "logit"),
    prior = rstanarm::normal(0, c(10, rep(5, dimension - 1L)), autoscale = FALSE),
    QR = FALSE, algorithm = "sampling", mean_PPD = FALSE,
    chains = chains, iter = nskip + retained, warmup = nskip, thin = 1,
    cores = 1, seed = seed, adapt_delta = 0.95, refresh = 0)

  # Extract post-warmup draws without shuffling chains. Missing coefficients
  # (for example, a predictor dropped by rstanarm) must cause an error.
  draws <- as.array(fit, pars = colnames(z))
  stopifnot(all(dim(draws) == c(retained, chains, dimension)),
            identical(dimnames(draws)[[3L]], colnames(z)))
  beta_draws <- matrix(draws, nrow = ndpost, ncol = dimension,
                       dimnames = list(NULL, colnames(z)))
  # R fills column by column: rows are chain 1, then chain 2, and so on.
  e_draws <- plogis(beta_draws %*% t(z)) #calculate propensity draws
  if (any(!is.finite(e_draws)) || any(e_draws <= 0 | e_draws >= 1))
    stop("Numerical propensity boundary reached; probabilities were not clipped.")

  # Stan's built-in checks exclude warmup. The divergence and maximum-tree-depth-hit counts below are for all chains in total;
  # E-BFMI measures energy exploration separately in each chain.
  info <- c(stan_seed = seed,
            stan_divergences = rstan::get_num_divergent(fit$stanfit),
            stan_max_treedepth_hits = rstan::get_num_max_treedepth(fit$stanfit),
            setNames(rstan::get_bfmi(fit$stanfit), paste0("stan_ebfmi_", seq_len(chains))))
  list(e = e_draws, beta = beta_draws, fit_info = info)
}
