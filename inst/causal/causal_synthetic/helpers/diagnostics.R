# posterior computes the diagnostics; this file prepares traces and flags concerns.
# No random numbers are drawn, and no draws or datasets are excluded.

diagnose_scalar <- function(values, chains, variable, component) {
  # Values contain complete chain blocks, so each matrix column is one chain.
  draws <- matrix(values, ncol = chains)
  finite <- all(is.finite(values))
  constant_chain <- FALSE
  for (chain in seq_len(ncol(draws))) {
    if (finite && length(unique(draws[, chain])) < 2L) constant_chain <- TRUE
  } 
  valid <- finite && !constant_chain && nrow(draws) >= 4L && chains >= 2L 

  # Unavailable diagnostics stay NA, including when a package call fails.
  metrics <- c(rhat_max = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
               ess_q025 = NA_real_, ess_q975 = NA_real_, mcse_mean = NA_real_,
               mcse_q025 = NA_real_, mcse_q975 = NA_real_)
  if (valid) {
    # Replace the values while keeping the metric names.
    metrics[] <- tryCatch( 
      suppressWarnings(c(
        posterior::rhat(draws),
        posterior::ess_bulk(draws),
        posterior::ess_tail(draws),
        posterior::ess_quantile(draws, probs = c(.025, .975), names = FALSE),
        posterior::mcse_mean(draws),
        posterior::mcse_quantile(draws, probs = c(.025, .975), names = FALSE)
      )),
      error = function(error) rep(NA_real_, length(metrics))
    )
  }

  # Express endpoint Monte Carlo errors as fractions of the interval width.
  endpoints <- c(NA_real_, NA_real_)
  if (finite) endpoints <- unname(quantile(values, c(.025, .975), type = 7))
  width <- endpoints[2L] - endpoints[1L]
  relative <- c(NA_real_, NA_real_)
  if (is.finite(width) && width > 0) {
    relative <- metrics[c("mcse_q025", "mcse_q975")] / width
  } 

  # Screening thresholds: missing diagnostics are flagged, never counted as a pass.
  rhat <- metrics[["rhat_max"]]
  effective_sizes <- metrics[c("ess_bulk", "ess_tail", "ess_q025", "ess_q975")]
  unavailable <- any(!is.finite(c(metrics, relative)))
  rhat_flag <- !is.finite(rhat) || rhat > 1.01 # we flag r hat values if diagnostic is unavailable or if rhat>1.01
  ess_flag <- any(!is.finite(effective_sizes)) || min(effective_sizes) < 100 * chains # we flag ESS if unavailable or the smallest of the four is below 100C, where C is number of chains
  endpoint_flag <- any(relative > .02, na.rm = TRUE) # flag if either endpoint's estimated MCSE exceeds 2% of interval width

  data.frame(
    variable = variable, component = component, as.list(metrics),
    q025 = endpoints[1L], q975 = endpoints[2L], interval_width = width,
    mcse_q025_relative_width = relative[1L], mcse_q975_relative_width = relative[2L],
    rhat_flag = rhat_flag, ess_flag = ess_flag,
    endpoint_mcse_flag = endpoint_flag, unavailable_flag = unavailable,
    diagnostic_flag = unavailable || rhat_flag || ess_flag || endpoint_flag,
    row.names = NULL, check.names = FALSE
  )
}

make_diagnostics <- function(x, a, y, m0, e, att, beta, chains) { # organizes the checks for one simulated dataset
  n <- length(a)
  number_draws <- length(att)
  stopifnot(
    length(chains) == 1L, is.finite(chains), chains >= 2L,
    chains == as.integer(chains), number_draws %% chains == 0L,
    length(y) == n, nrow(x) == n, all(a %in% c(0L, 1L)),
    any(a == 0L), any(a == 1L), # requires at least one treated and at least one control individual
    identical(dim(m0), c(number_draws, n)), identical(dim(e), c(number_draws, n)),
    nrow(beta) == number_draws, ncol(beta) == ncol(x) + 1L
  )

  # 1. Assemble the monitored quantities: one row per draw, one column per quantity.
  treated <- a == 1L
  control <- !treated
  number_treated <- sum(treated)
  # Outcome-model diagnostics monitor the mean m0 prediction among treated
  # observations and among controls, rather than every individual prediction.
  # The ATT before and after Bayesian-bootstrap weighting is also monitored below.
  labels <- c("m0_mean_treated", "m0_mean_control", "e_mean_treated", "e_mean_control",
              "max_control_odds", "att_without_bb") # six quantities that will be constructed from each pair of nuisance-model draws
  values <- matrix(NA_real_, number_draws, length(labels), dimnames = list(NULL, labels))
  for (draw in seq_len(number_draws)) {
    residual <- y - m0[draw, ]
    odds <- e[draw, control] / (1 - e[draw, control])
    # Equal weights track model draws before adding Bayesian-bootstrap randomness.
    before_bb <- (sum(residual[treated]) - sum(odds * residual[control])) / number_treated
    values[draw, ] <- c(mean(m0[draw, treated]), mean(m0[draw, control]),
                       mean(e[draw, treated]), mean(e[draw, control]),
                       max(odds), before_bb)
  }
  colnames(beta) <- paste0("alpha_", 0:(ncol(beta) - 1L)) # name the propensity coefficient columns
  values <- cbind(values, corrected_att = att, beta)

  # 2. Label chain blocks; paired BART and propensity chains remain separate samplers.
  per_chain <- number_draws %/% chains # number of draws per chain
  trace <- data.frame(chain = rep(seq_len(chains), each = per_chain),
                      iteration = rep(seq_len(per_chain), times = chains),
                      values, check.names = FALSE)

  # 3. Group the quantities by model, then diagnose each column in order.
  variables <- colnames(values) 
  components <- variables 
  components[startsWith(variables, "m0_")] <- "outcome" 
  components[startsWith(variables, "e_") | variables == "max_control_odds"] <- "propensity"
  components[startsWith(variables, "alpha_")] <- "logistic_parameters"
  diagnostic_rows <- vector("list", ncol(values))
  for (column in seq_len(ncol(values))) {
    diagnostic_rows[[column]] <- diagnose_scalar(
      values[, column], chains, variables[column], components[column]
    )
  }
  list(trace = trace, table = do.call(rbind, diagnostic_rows))
}
