## this script comes from the unbiasedmcmc repository on github
## https://raw.githubusercontent.com/pierrejacob/unbiasedmcmc/master/inst/reproduceepidemiology/plummer.cut.run.R

## this script is way too complicated for what it does; couplings are not really used in the end; to simplify

# load packages
library(unbiasedmcmc)
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

sample_module1 <- function(nsamples){
  theta1s <- matrix(nrow = nsamples, ncol = J)
  for (j in 1:J){
    theta1s[,j] <- rbeta(nsamples, shape1 = posterior_phi_alpha[j], shape2 = posterior_phi_beta[j])
  }
  return(theta1s)
}

get_kernels <- function(theta1, Sigma_proposal, init_mean, init_Sigma){
  target <- function(x) epidemiology_loglikC(theta1, x, ncases, Npop_normalized) + dprior2(x, hyper2)
  ##
  Sigma_proposal_chol <- chol(Sigma_proposal)
  Sigma_proposal_chol_inv <- solve(Sigma_proposal_chol)

  # Markov kernel of the chain
  single_kernel <- function(state){
    chain_state <- state$chain_state
    current_pdf <- state$current_pdf
    proposal_value <- chain_state + fast_rmvnorm_chol(1, rep(0, dimension), Sigma_proposal_chol)[1,]
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }

  # Markov kernel of the coupled chain
  coupled_kernel <- function(state1, state2){
    chain_state1 <- state1$chain_state
    chain_state2 <- state2$chain_state
    current_pdf1 <- state1$current_pdf
    current_pdf2 <- state2$current_pdf
    proposal_value <- rmvnorm_reflectionmax(chain_state1, chain_state2, Sigma_proposal_chol, Sigma_proposal_chol_inv)
    proposal1 <- proposal_value$xy[,1]
    proposal2 <- proposal_value$xy[,2]
    identical_proposal <- proposal_value$identical
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- proposal_pdf1
    if (!identical_proposal){
      proposal_pdf2 <- target(proposal2)
    }
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    identical_ <- ((proposal_value$identical) && (accept1) && (accept2))
    return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
                state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
                identical = identical_))
  }
  rinit <- function(){
    chain_state <- fast_rmvnorm_chol(1, mean = init_mean, init_Sigma)[1,]
    current_pdf <- target(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  return(list(target = target, coupled_kernel = coupled_kernel, single_kernel = single_kernel, rinit = rinit))
}

# (modify nsamples if desired)
# first sample theta1's and compute posterior mean under first model
nsamples <- 1000
theta1s <- sample_module1(nsamples)
theta1hat <- colMeans(theta1s)
## then try to perform inference on theta2 given theta1hat
dimension <- 2
Sigma_proposal <- diag(0.1, dimension, dimension)
init_mean <- rep(0, dimension)
init_Sigma <- diag(1, dimension, dimension)
## first run standard MCMC
pb <- get_kernels(theta1hat, Sigma_proposal, init_mean, init_Sigma)
niterations <- 5e3
chain <- matrix(0, nrow = niterations, ncol = dimension)
state <- pb$rinit()
for (iteration in 1:niterations){
  state <- pb$single_kernel(state)
  chain[iteration,] <- state$chain_state
}
# matplot(chain, type = "l")
# then refine the initial distribution and the estimate of posterior covariance matrix
chain_postburn <- chain[2e3:niterations,]
init_mean <- colMeans(chain_postburn)
init_Sigma <- cov(chain_postburn)
Sigma_proposal <- init_Sigma
pb <- get_kernels(theta1hat, Sigma_proposal, init_mean, init_Sigma)
chain <- matrix(0, nrow = niterations, ncol = dimension)
state <- pb$rinit()
for (iteration in 1:niterations){
  state <- pb$single_kernel(state)
  chain[iteration,] <- state$chain_state
}
# matplot(chain, type = "l")
chain_postburn <- chain[2e3:niterations,]
init_mean <- colMeans(chain_postburn)
init_Sigma <- cov(chain_postburn)
Sigma_proposal <- init_Sigma
# the chain seems to mix OK

## Now let theta1 vary too and see how
## meeting times behave with the tuning selected as above
theta1s <- sample_module1(nsamples)
meetingtimes_1 <-  foreach(irep = 1:nsamples) %dorng% {
  pb <- get_kernels(theta1s[irep,], Sigma_proposal, init_mean, init_Sigma)
  sample_meetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)
}
meetingtime <- sapply(meetingtimes_1, function(x) x$meetingtime)
summary(meetingtime)
## the meeting times are quite short
k <- 100
m <- 1000

# now we will re-estimate the target covariance matrix using unbiased MCMC
c_chains_1 <-  foreach(irep = 1:nsamples) %dorng% {
  pb <- get_kernels(theta1s[irep,], Sigma_proposal, init_mean, init_Sigma)
  sample_coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
#
max(sapply(c_chains_1, function(x) x$meetingtime))

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_1[[irep]], h = function(x) x, k = k, m = m)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_1[[irep]], h = function(x) x^2, k = k, m = m)
}

cross_estimator <- foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_1[[irep]], h = function(x) x[1]*x[2], k = k, m = m)
}

post_mean <- rowMeans(sapply(mean_estimators, function(x) x))
post_var <- rowMeans(sapply(square_estimators, function(x) x)) - post_mean^2
post_cross <- mean(sapply(cross_estimator, function(x) x)) - prod(post_mean)
Sigma_proposal <- diag(post_var)
Sigma_proposal[1,2] <- Sigma_proposal[2,1] <- post_cross
init_Sigma <- Sigma_proposal
## we could check that the estimated covariance is positive semi-definite, e.g. by running
solve(Sigma_proposal)
## with the above tuning parameters, we are now ready to produce the final results

# ## Using the new covariance matrix estimate we draw new meeting times
nsamples <- 10000
theta1s <- sample_module1(nsamples)

### exact cut distribution from parallel MCMC runs
## Modify niterations
niterations <- 1000
theta2s <- foreach(itheta = 1:nrow(theta1s), .combine = rbind) %dorng% {
  theta1 <- theta1s[itheta,]
  pb <-  get_kernels(theta1, Sigma_proposal, init_mean, init_Sigma)
  chain <- pb$rinit()
  for (iter in 1:niterations){
    chain <- pb$single_kernel(chain)
  }
  chain$chain_state
}
save(theta2s, file = "inst/epidemio/epidemio_cut_mcmc_samples.RData")
# load(file = "plummer.mcmc.RData")
