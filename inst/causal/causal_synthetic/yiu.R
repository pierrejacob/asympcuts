# [Y1] Our observed data replace ACIC; [Y2] binary outcome BART;
# [Y3] logistic propensity posterior; [Y4] authors' correction unchanged;
# [Y5] diagnostics; [Y6] raw ATT draws (Python summarizes once).
yiu_att_fit <- function(data, seed, package_dir, ndpost, nskip, chains, ntree) {
  stopifnot(as.character(utils::packageVersion("BART")) == "2.9.10",
            requireNamespace("posterior", quietly = TRUE),
            requireNamespace("rstanarm", quietly = TRUE), .Platform$OS.type == "unix") # mc.pbart() rejects non-Unix systems
  controls <- c(ndpost, nskip, chains, ntree, seed)
  stopifnot(all(is.finite(controls)), all(controls == floor(controls)),
            nskip >= 0, chains >= 2, ntree > 0, seed > 0,
            seed <= .Machine$integer.max - 2, ndpost %% chains == 0,
            ndpost / chains >= 20, chains <= parallel::detectCores())
  x <- as.matrix(data$x)
  a <- data$a
  y <- data$y
  n <- nrow(x)
  stopifnot(is.numeric(x), ncol(x) == 9, all(is.finite(x)), length(a) == n, length(y) == n,
            all(a %in% c(0, 1)), all(y %in% c(0, 1)),
            length(unique(a)) == 2, length(unique(y)) == 2)
  a <- as.integer(a)
  y <- as.integer(y)

  original <- new.env()  # Keep the authors' functions in their own environment, which is hereby created.
  original_dir <- file.path(package_dir, "original_yiu", "ACIC2016")
  sys.source(file.path(original_dir, "util_fun.R"), envir = original)
  sys.source(file.path(original_dir, "1step_ATT_fun.R"), envir = original)
  source(file.path(package_dir, "helpers", "logistic_posterior.R"), local = TRUE)
  source(file.path(package_dir, "helpers", "diagnostics.R"), local = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection") # specify the random-number algorithms

  # [Y2] Pool both treatment arms; prediction rows first set A=1, then A=0.
  outcome <- BART::mc.pbart(
    x.train = cbind(x, treatment = a), y.train = y,
    x.test = rbind(cbind(x, treatment = 1L), cbind(x, treatment = 0L)),
    ntree = ntree, ndpost = ndpost, nskip = nskip, mc.cores = chains, 
    keepevery = 1L, nice = 19L, seed = seed, printevery = ndpost + nskip + 1L)
  stopifnot(all(dim(outcome$prob.test) == c(ndpost, 2 * n))) #check the returned dimensions: one row per retained draw and 2n columns, i.e. one per prediction input
  m1 <- outcome$prob.test[, seq_len(n), drop = FALSE] # treated prediction
  m0 <- outcome$prob.test[, n + seq_len(n), drop = FALSE] # untreated prediction
  stopifnot(all(is.finite(m0)), all(is.finite(m1)),
            all(m0 >= 0 & m0 <= 1), all(m1 >= 0 & m1 <= 1))
  rm(outcome) 
  invisible(gc(verbose = FALSE)) #release unused memory

  # [Y3] Outcome Y is not an input to the propensity posterior.
  propensity <- fit_logistic_posterior(cbind(intercept = 1, x), a,
                                        ndpost, nskip, chains, seed + 1L)
  e <- propensity$e
  # [Y4] Every call draws fresh Dirichlet weights inside original onestep().
  set.seed(seed + 2L)
  att <- numeric(ndpost)
  for (draw in seq_len(ndpost)) {
    att[draw] <- original$onestep(fdraw0 = m0[draw, ], fdraw1 = m1[draw, ],
                                pi = e[draw, ], r = a, y = y, n = n)[2L] # calls on the onestep function from authors' original implementation to obtain the att
  }
  stopifnot(all(is.finite(att)))
  # [Y5] Diagnostics may flag concerns but never remove draws or refit a model.
  extra <- make_diagnostics(x, a, y, m0, e, att, propensity$beta, chains)
  info <- c(seed = seed, propensity$fit_info) #retaining Stan diagnostics
  list(draws = data.frame(att = att), trace = extra$trace, diagnostics = extra$table,
       fit_info = data.frame(item = names(info), value = unname(info))) 
}

# This script runs in a fresh R process.
if (sys.nframe() == 0L) { #should only run when we execute the script directly
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 7L)
    stop("Use: Rscript yiu.R OBSERVED OUTPUT SEED NDPOST NSKIP CHAINS NTREE")
  script <- grep("^--file=", commandArgs(), value = TRUE)[1L] #find the script's own path
  package_dir <- dirname(normalizePath(sub("^--file=", "", script)))
  observed <- read.csv(args[1L], check.names = FALSE)
  stopifnot(ncol(observed) == 11L, identical(tail(names(observed), 2L), c("A", "Y")))
  fit <- yiu_att_fit(list(x = as.matrix(observed[, 1:9]), a = observed$A, y = observed$Y),
                     as.numeric(args[3L]), package_dir, as.integer(args[4L]),
                     as.integer(args[5L]), as.integer(args[6L]), as.integer(args[7L]))
  dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
  # [Y6] Each quantity is stored once; Python computes intervals from yiu_draws.csv.
  files <- c(draws = "yiu_draws", trace = "trace", diagnostics = "diagnostics", fit_info = "fit_info")
  for (name in names(files))
    write.csv(fit[[name]], file.path(args[2L], paste0(files[[name]], ".csv")), row.names = FALSE)
}
