
## using the generic RWMRTH implementation
target <- list(dimension = J, density = function(theta1){
  return(sapply(1:dim(theta1)[1], function(ichain) logposterior1(theta1[ichain,])))
})

## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(rep(0.01, J))
tuning_parameters$niterations <- 100000
tuning_parameters$nchains <- 4
tuning_parameters$adaptation <- 1000
tuning_parameters$rinit <- function(nchains){
  chains <- matrix(runif(nchains*target$dimension, 0, 1), ncol = target$dimension)
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
dim(mcmc_samples_stage1)



# library(gridExtra)
# df_stage1 <- data.frame(method = "mcmc", theta = mcmc_samples_stage1)
# g1 <- ggplot(df_stage1, aes(x = theta.1, fill = method)) + geom_density(alpha = 0.5)
# ## add beta posterior density curve
# g1 <- g1 + stat_function(fun = dbeta, args = list(shape1 = posterior_phi_alpha[1],
#                                      shape2 = posterior_phi_beta[1]))
# g2 <- ggplot(df_stage1, aes(x = theta.2, fill = method)) + geom_density(alpha = 0.5)
# g2 <- g2 + stat_function(fun = dbeta, args = list(shape1 = posterior_phi_alpha[2],
#                                      shape2 = posterior_phi_beta[2]))
# g3 <- ggplot(df_stage1, aes(x = theta.3, fill = method)) + geom_density(alpha = 0.5)
# g3 <- g3 + stat_function(fun = dbeta, args = list(shape1 = posterior_phi_alpha[3],
#                                      shape2 = posterior_phi_beta[3]))
# grid.arrange(g1, g2, g3,  ncol = 3)

##
N1 <- 5000
selected_indices <- round(seq(1, nrow(mcmc_samples_stage1), length.out = N1))

rwmrth_res_2 <- mclapply(selected_indices, function(i){
  some_theta1 <- mcmc_samples_stage1[i,]
  ## using the generic RWMRTH implementation
  target2given1 <- list(dimension = 2, density = function(theta2){
    return(sapply(1:dim(theta2)[1], function(ichain) logposterior2(some_theta1, theta2[ichain,])))
  })
  ## tuning parameters for rwmrth
  tuning_parameters <- list()
  tuning_parameters$cov_proposal <- diag(c(1, 1))
  ## length of each chain
  tuning_parameters$niterations <- 5000
  ## number of chains
  tuning_parameters$nchains <- 2
  ## adapt covariance every 1000 iterations
  tuning_parameters$adaptation <- 1000
  tuning_parameters$rinit <- function(nchains){
    return(matrix(rnorm(nchains*target2given1$dimension, 0, 1), ncol = target2given1$dimension))
  }
  ## initialization
  current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
  ## MCMC run
  rwmrth_res_2 <- rwmrth(target2given1, tuning_parameters)
  ## default choice of burnin
  ## trace plot
  burnin <- tuning_parameters$niterations/2
  # matplot(rwmrth_res_2$chains[[1]][burnin:tuning_parameters$niterations,], type = 'l')
  # take last 2 samples from each chain
  mcmc_samples_stage2 <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
    rwmrth_res_2$chains[[ichain]][(tuning_parameters$niterations-1):tuning_parameters$niterations,]
  }
  mcmc_samples_stage2
}, mc.cores = 10)

theta1s <- mcmc_samples_stage1[selected_indices,]
theta2s <- do.call('rbind', rwmrth_res_2 ) %>% as.matrix()
save(theta1s, theta2s, file = "inst/epidemio/epidemio_cut_mcmc_samples.RData")

# df_mcmc <- do.call('rbind', rwmrth_res_2) %>% as.data.frame()
# df_mcmc$method <- 'cut Bayesian'
# dim(df_mcmc)
# head(df_mcmc)
#
# ## plot 95% ellipse from V1,V2
# library(ggplot2)
# g_cut <- ggplot(df_mcmc, aes(x = V1, y = V2)) +
#   stat_ellipse(level = 0.95, type = 'norm', size = 1) +
#   xlab(expression(theta[2.1])) + ylab(expression(theta[2.2]))
# g_cut

