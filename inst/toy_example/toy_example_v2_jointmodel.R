#### Data ####
load("inst/toy_example/toy_example_data.RData")
if (SCENARIO == 1){
  attach(data_scenario1)
}
if(SCENARIO == 2){
  attach(data_scenario2)
}
#### Log-posterior density function and point estimation ####
log_posterior <- function(theta){
  return(ll1(z, theta[1]) + ll2(theta[1], theta[2], x, y) + lp1(theta[1]) + lp2(theta[2]))
}
## joint model's posterior for theta_1 and theta_2
target <- list(dimension = 2, density = function(theta){
  return(sapply(1:dim(theta)[1], function(ichain) log_posterior(theta[ichain,])))
})

##### joint MLE #####
# Function to optimize
optim_ <- tryCatch(optim(par = c(0,0),
                               fn = function(th) - log_posterior(th),
                               method = 'BFGS'),
                         error = function(e)  return(list(error = TRUE)))
jointmodel_map <- optim_$par
jointmodel_map[1]
# nearly identical to
mean((z+y)/2) - jointmodel_map[2]/2 * mean(x)
# or, purely in terms of the data
{sum(z+y) - sum(x) * sum(x*y)/sum(x^2)}/{2*n*(1 - (sum(x)^2)/(2*n*sum(x^2)))}
# limiting value of 1st component of joint model MLE:
meanz = 0
meanx = 1
varx = 1
meany = meanx + 1
meanx2 = meanx^2 + varx
meanxy = meanx2 + meanx
0.5 * (meanz + meany - meanx * meanxy / meanx2) / (1 - meanx^2 / (2 * meanx2))

# and for the second parameter
jointmodel_map[2]
sum(x*y)/sum(x^2) - sum(x) / (sum(x^2)) * jointmodel_map[1]

#### Joint model Bayesian inference via RWMRTH ####

## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(c(1,1))
tuning_parameters$niterations <- 10000
tuning_parameters$nchains <- 4
tuning_parameters$adaptation <- 1000
tuning_parameters$rinit <- function(nchains){
  matrix(rnorm(nchains*target$dimension, 0, 1), ncol = target$dimension)
}
## initialization
current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
## MCMC run
rwmrth_res_joint <- rwmrth(target, tuning_parameters)
## default choice of burnin
burnin <- tuning_parameters$niterations/5
## trace plot of first column of each matrix in rwmrth_res_joint$chains
# plot(rwmrth_res_joint$chains[[1]][,1], type ='l', ylim = c(-2, +2))
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][,1])
# # zoom in
# plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,1], type ='l', ylim = c(-.1, .5), ylab = expression(theta[1]))
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,1])
# ## same thing for second component
# plot(rwmrth_res_joint$chains[[1]][,2], type ='l', ylim = c(0, +2))
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][,2])
# # zoom in
# plot(rwmrth_res_joint$chains[[1]][burnin:tuning_parameters$niterations,2], type ='l', ylim = c(1.2, 1.6), ylab = expression(theta[2]))
# for (ichain in 2:tuning_parameters$nchains) lines(rwmrth_res_joint$chains[[ichain]][burnin:tuning_parameters$niterations,2])
## scatter plot of first two columns, after burn-in
mcmc_samples <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
  rwmrth_res_joint$chains[[ichain]][(burnin+1):tuning_parameters$niterations,]
}
# plot(mcmc_samples[,1], mcmc_samples[,2], pch = 16, cex = 0.5, xlab = expression(theta[1]), ylab = expression(theta[2]))
# ## add joint model MLE
# points(jointmodel_map[1], jointmodel_map[2], col ='red', pch = 19, cex = 3)
#
# colMeans(mcmc_samples)
# cov(mcmc_samples)

#### Joint model Laplace approximation ####

# Compute Hessian at joint MLE
H <- -hessian(func = log_posterior, jointmodel_map)
cov_laplace <- solve(H)
samples_laplace <- fast_rmvnorm(1e5, jointmodel_map, cov_laplace)

#### Plot Ellipses for the joint modelling approach ####
## plot 95% credible region associated with mcmc_samples using stat_ellipse
# df_samples <- data.frame(theta1 = mcmc_samples[,1], theta2 = mcmc_samples[,2])
# ggplot(df_samples, aes(x = theta1, y = theta2)) + theme_bw() + xlim(0.2,0.5) + ylim(1.2, 1.5) +
#   stat_ellipse(level = 0.95, color = 'blue', type = 'norm', linewidth = 2) +
#   stat_ellipse(level = 0.95, type = 'norm', mapping = aes(x = theta1, y = theta2),
#                data = data.frame(theta1 = samples_laplace[,1], theta2 = samples_laplace[,2]),
#                inherit.aes = FALSE,
#                color = 'red', linetype = 2, linewidth = 2) +
#   geom_point(data = data.frame(theta1 = jointmodel_map[1], theta2 = jointmodel_map[2]), aes(x = theta1, y = theta2),
#              color ='black', pch = 19, size = 3) +
#   xlab(expression(theta[1])) + ylab(expression(theta[2])) +
#   ggtitle('95% credible regions: Standard Posterior (blue) vs Laplace approximation (red)')

save(mcmc_samples, jointmodel_map, cov_laplace, file = paste0('inst/toy_example/toy_example_v2_jointmodel_scenario', SCENARIO, '.RData'))
