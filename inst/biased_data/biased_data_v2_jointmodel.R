attach(biased_dataset)


#### Standard Bayesian inference ####
#### Log-posterior density function and point estimation ####
## write log-likelihood function of the joint model for (theta_1,theta_2)
log_lik_joint <- function(theta){
  part1 <- sum(dnorm(x1, mean = theta[1], sd = 1, log = TRUE))
  part2 <- sum(dnorm(x2, mean = theta[1] + theta[2], sd = 1, log = TRUE))
  return(part1 + part2)
}
log_prior_density <- function(theta){
  prior1 <- dnorm(theta[1], mean = 0, sd = sigma_1_prior, log = TRUE)
  prior2 <- dnorm(theta[2], mean = 0, sd = sigma_2_prior, log = TRUE)
  return(prior1 + prior2)
}

## joint model's posterior for theta_1 and theta_2
target_density <- function(theta) log_lik_joint(theta) + log_prior_density(theta)
## target object for MCMC sampling (where theta is a matrix with each row a different chain)
target <- list(dimension = 2, density = function(theta){
  return(sapply(1:dim(theta)[1], function(ichain) target_density(theta[ichain,])))
})

#### MCMC for joint model Bayesian inference ####

## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(c(1,1))
tuning_parameters$niterations <- 5000
tuning_parameters$nchains <- 4
tuning_parameters$adaptation <- 1000
tuning_parameters$rinit <- function(nchains){
  matrix(rnorm(nchains*target$dimension, 0, 1), ncol = target$dimension)
}
## initialize chains
current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
## run rwmrth
rwmrth_res_joint <- rwmrth(target, tuning_parameters)
burnin <- tuning_parameters$niterations/5
## trace plot of first column of each matrix in rwmrth_res_joint$chains
# plot(rwmrth_res_joint$chains[[1]][,1], type ='l')
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][,1])
# # zoom in
# plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,1], type ='l', ylab = expression(theta[1]))
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,1])
# ## same thing for second component
# plot(rwmrth_res_joint$chains[[1]][,2], type ='l')
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][,2])
# # zoom in
# plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,2], type ='l', ylab = expression(theta[2]))
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,2])
## scatter plot of first two columns, after burn-in
mcmc_samples <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
  rwmrth_res_joint$chains[[ichain]][(burnin+1):tuning_parameters$niterations,]
}

#### Joint model Laplace approximation ####
optim_ <- tryCatch(optim(par = c(0,0),
                         fn = function(x) -target_density(x),
                         method = 'BFGS'),
                   error = function(e)  return(list(error = TRUE)))
jointmodel_map <- optim_$par
# plot(mcmc_samples[,1], mcmc_samples[,2], pch = 16, cex = 0.5, xlab = expression(theta[1]), ylab = expression(theta[2]))
# ## add joint model MAP
# points(jointmodel_map[1], jointmodel_map[2], col ='red', pch = 19, cex = 3)

H <- -hessian(func = target_density, jointmodel_map)
cov_laplace <- solve(H)

# samples_laplace <- fast_rmvnorm(1e5, jointmodel_map, cov_laplace)

#### Plot Ellipses for the joint modelling approach ####
## plot 95% credible region associated with mcmc_samples using stat_ellipse
# df_samples <- data.frame(theta1 = mcmc_samples[,1], theta2 = mcmc_samples[,2])
# ggplot(df_samples, aes(x = theta1, y = theta2)) + theme_bw() + #xlim(0.2,0.5) + ylim(1.2, 1.5) +
#   stat_ellipse(level = 0.95, color = 'blue', type = 'norm', linewidth = 2) +
#   stat_ellipse(level = 0.95, type = 'norm', mapping = aes(x = theta1, y = theta2),
#                data = data.frame(theta1 = samples_laplace[,1], theta2 = samples_laplace[,2]),
#                inherit.aes = FALSE,
#                color = 'red', linetype = 2, linewidth = 2) +
#   xlab(expression(theta[1])) + ylab(expression(theta[2])) +
#   ggtitle('95% credible regions: Standard Posterior (blue) vs Laplace approximation (red)')

save(mcmc_samples, jointmodel_map, cov_laplace, file = paste0("inst/biased_data/biased_data_v2_jointmodel_n1_", n1, ".RData"))
detach(biased_dataset)

