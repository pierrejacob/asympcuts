#'@rdname rwmrth
#'@title Random walk Metropolis-Rosenbluth-Teller-Hastings algorithm
#'@description with multiple chains, adaptive proposal
#' @param tuning_parameters must be a list containing:
#' - niterations: number of iterations
#' - nchains: number of chains
#' - cov_proposal: initial proposal covariance matrix
#' - rinit: function for initializing the chains; it must take as input
#' the number of chains and return a matrix with nchains rows and target$dimension columns
#' - adaptation: adaptation frequency (if 0 no adaptation)
#' @param target target distribution; it must be a list containing:
#' - dimension: dimension of the target
#' - density: function for evaluating the log target density; it must take as input
#' a matrix with nchains rows and target$dimension columns and return a vector of length nchains
#' @return a list containing:
#' - chains: a list of length nchains; each element is a matrix with niterations rows and target$dimension columns
#' - current_dtarget: vector of length nchains with the log target density values at the final states of the chains
#' - dtarget_history: matrix with niterations rows and nchains columns with the log target density values at each iteration
#' - naccepts: total number of accepted proposals across all chains
#' - cov_proposal: final proposal covariance matrix
#'@export
rwmrth <- function(target, tuning_parameters){
  niterations <- tuning_parameters$niterations
  nchains <- tuning_parameters$nchains
  cov_proposal <- tuning_parameters$cov_proposal
  # store whole chains
  chains <- rep(list(matrix(nrow = niterations, ncol = target$dimension)), nchains)
  # current states of the chains
  current_chains <- matrix(nrow = nchains, ncol = target$dimension)
  # initialization of the chains
  current_chains <- tuning_parameters$rinit(nchains)
  for (ichain in 1:nchains){
    chains[[ichain]][1,] <- current_chains[ichain,]
  }
  # log target density values associated with the current states of the chains
  current_dtarget <- target$density(current_chains)
  dtarget_history <- matrix(nrow = niterations, ncol = nchains)
  dtarget_history[1,] <- current_dtarget
  #
  naccepts <- 0
  # run the chains
  for (iteration in 2:niterations){
    if (iteration > 50 && tuning_parameters$adaptation > 0  && (iteration %% tuning_parameters$adaptation) == 0){
      # adapt the proposal covariance matrix based on the last < 50,000 samples of all chains
      mcmc_samples <- foreach(ichain = 1:nchains, .combine = rbind) %do% {
        matrix(chains[[ichain]][max(iteration-tuning_parameters$adaptation, iteration - 50000):(iteration-1),], ncol = target$dimension)
      }
      cov_proposal <- cov(mcmc_samples) / target$dimension
    }
    # proposals
    proposals <- current_chains + fast_rmvnorm(nchains, rep(0, target$dimension), cov_proposal)
    # proposals' target density
    proposal_dtarget <- target$density(proposals)
    # log Metropolis Hastings ratio
    acceptance_ratios <- (proposal_dtarget - current_dtarget)
    # uniforms for the acceptance decisions
    uniforms <- runif(n = nchains)
    # acceptance decisions
    accepts <- (log(uniforms) < acceptance_ratios)
    naccepts <- naccepts + sum(accepts)
    # make the appropriate replacements
    current_chains[accepts,] <- proposals[accepts,]
    if (is.null(dim(current_chains))) current_chains <- matrix(current_chains, ncol = target$dimension)
    current_dtarget[accepts] <- proposal_dtarget[accepts]
    # book keeping
    for (ichain in 1:nchains){
      chains[[ichain]][iteration,] <- current_chains[ichain,]
    }
    dtarget_history[iteration,] <- current_dtarget
  }
  cat("average acceptance:", naccepts / (niterations*nchains) * 100, "%\n")
  return(list(chains = chains, current_dtarget = current_dtarget, dtarget_history = dtarget_history,
              naccepts = naccepts, cov_proposal = cov_proposal))
}
