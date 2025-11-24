### runs RWMRTH for joint model's posterior distribution in Plummer's epidemiological model

##

## joint model's posterior for theta_1 and theta_2
target <- list(dimension = 15, density = function(theta){
  return(sapply(1:dim(theta)[1], function(ichain) logposterior1(theta[ichain,1:13]) + logposterior2(theta[ichain,1:13], theta[ichain,14:15])))
})



## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(rep(0.1, 15))
tuning_parameters$niterations <- 50000
tuning_parameters$nchains <- 4
tuning_parameters$adaptation <- 1000
tuning_parameters$rinit <- function(nchains){
  chains <- matrix(rnorm(nchains*target$dimension, 0, 1), ncol = target$dimension)
  for (j in 1:13){
    chains[,j] <- rbeta(nchains, shape1 = posterior_phi_alpha[j], shape2 = posterior_phi_beta[j])
  }
  return(chains)
}
chains_ <- tuning_parameters$rinit(2)
target$density(chains_)
## initialization
current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
## MCMC run
rwmrth_res_joint <- rwmrth(target, tuning_parameters)
## default choice of burnin
burnin <- tuning_parameters$niterations/5

# ## trace plot of first column of each matrix in rwmrth_res_joint$chains
# par(mfrow = c(2,3))
# for (i in 1:6){
#   plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,i], type ='l')
#   for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,i])
# }

##
# par(mfrow = c(3,1))
# plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,13], type ='l')
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,13])
# for (i in 14:15){
#   plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,i], type ='l')
#   for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,i])
# }

joint_mcmc_samples <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
  rwmrth_res_joint$chains[[ichain]][(burnin+1):tuning_parameters$niterations,]
}

# par(mfrow = c(1,1))
# matplot(joint_mcmc_samples[,1:4], type = 'l')
# matplot(joint_mcmc_samples[,14], type = 'l')
# matplot(joint_mcmc_samples[,15], type = 'l')

save(tuning_parameters, burnin, joint_mcmc_samples, file = 'inst/epidemio/epidemio_joint_mcmc_samples.RData')


