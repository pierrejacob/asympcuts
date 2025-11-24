# load data
attach(biased_dataset)
#### Two-stage inference ####

target_density1 <- function(theta1) sum(dnorm(x1, mean = theta1, sd = 1, log = TRUE)) + dnorm(theta1, mean = 0, sd = sigma_1_prior, log = TRUE)

## first target
target1 <- list(dimension = 1, density = function(theta1){
  return(sapply(1:dim(theta1)[1], function(ichain) target_density1(theta1[ichain,1])))
})

# get second stage target given theta1
get_target2 <- function(theta1){
  target_density2_ <- function(theta2) sum(dnorm(x2, mean = theta1 + theta2, sd = 1, log = TRUE)) + dnorm(theta2, mean = 0, sd = sigma_2_prior, log = TRUE)
  ## joint model's posterior for theta_1 and theta_2
  return(list(dimension = 1, density = function(theta2){
    return(sapply(1:dim(theta2)[1], function(ichain)  target_density2_(theta2[ichain,1])))
  }))
}

##### 2SMAP #####
laplace_cut_res <- laplace_cut(
  logposterior1 = target_density1,
  logposterior2 = function(theta1, theta2) {
    sum(dnorm(x2, mean = theta1 + theta2, sd = 1, log = TRUE)) + dnorm(theta2, mean = 0, sd = sigma_2_prior, log = TRUE)
  },
  init_theta1 = 0,
  init_theta2 = 0
)

laplace_cut_res$thetahat[1]
# nearly identical to
mean(x1)
# and second stage MLE
laplace_cut_res$thetahat[2]
sum(x2)/n2 - laplace_cut_res$thetahat[1]
# print summary:
cat("Two-stage MAP:\n", laplace_cut_res$thetahat, "\n")

#### Cut Bayesian inference via RWMRTH ####
##### Stage 1 #####
## tuning parameters for rwmrth
tuning_parameters <- list()
tuning_parameters$cov_proposal <- diag(1)
tuning_parameters$niterations <- 5000
tuning_parameters$nchains <- 4
tuning_parameters$adaptation <- 1000
tuning_parameters$rinit <- function(nchains){
  matrix(rnorm(nchains*target1$dimension, 0, 1), ncol = target1$dimension)
}
current_chains <- tuning_parameters$rinit(tuning_parameters$nchains)
rwmrth_stage1_res <- rwmrth(target1, tuning_parameters)
burnin <- tuning_parameters$niterations/5
#
all_samples1 <- foreach(ichain = 1:tuning_parameters$nchains, .combine = rbind) %do% {
  rwmrth_stage1_res$chains[[ichain]][(burnin+1):tuning_parameters$niterations,,drop=F]
}

# compare to Laplace approximation
# Compute Hessian at 1st-stage MLE
var(all_samples1[,1])
laplace_cut_res$asympvar[1,1]
## very accurate

##### Stage 2 #####
# select nsamples equispaced draws in all_samples1
nsamples1 <- 5000
# stop if nsamples1 < nrow(all_samples1)
if (nsamples1 > nrow(all_samples1)) {
  stop("nsamples1 is larger than the number of available samples from stage 1")
}

selected_indices <- round(seq(1, nrow(all_samples1), length.out = nsamples1))
theta1_samples_for_stage2 <- all_samples1[selected_indices, , drop=F]
## register parallel backend
library(doParallel)
library(doRNG)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("biased_dataset"))
clusterEvalQ(cl, {
  library(asympcuts)
})


# foreach loop, loading functions and package in memory
theta2_samples <- foreach(isample1 = 1:nsamples1, .combine = rbind) %dorng% {
  attach(biased_dataset)
  theta1_sample <- theta1_samples_for_stage2[isample1, 1]
  target2 <- get_target2(theta1_sample)
  # tuning parameters for rwmrth
  tuning_parameters_stage2 <- list()
  tuning_parameters_stage2$cov_proposal <- diag(1)
  tuning_parameters_stage2$niterations <- 300
  tuning_parameters_stage2$nchains <- 1
  tuning_parameters_stage2$adaptation <- 100
  tuning_parameters_stage2$rinit <- function(nchains){
    matrix(rnorm(nchains*target2$dimension, 0, 1), ncol = target2$dimension)
  }
  rwmrth_stage2_res <- rwmrth(target2, tuning_parameters_stage2)
  # plot(rwmrth_stage2_res$chains[[1]][,1], type ='l')
  detach(biased_dataset)
  # collect last sample
  tail(rwmrth_stage2_res$chains[[1]], n = 1)[1]
}
# end cluster
stopCluster(cl)

cut_samples <- cbind(theta1_samples_for_stage2, as.numeric(theta2_samples))


plot(cut_samples[,1], cut_samples[,2], pch = 16, cex = 0.5, xlab = expression(theta[1]), ylab = expression(theta[2]))
## add 2S MAP
points(laplace_cut_res$thetahat[1], laplace_cut_res$thetahat[2], col ='red', pch = 19, cex = 3)

#### Cut Laplace approximation ####

## Two stage MAP obtained above
twostage_map <- laplace_cut_res$thetahat
colMeans(cut_samples)
twostage_map

##
cov(cut_samples)
laplace_cut_res$asympvar

cut_laplace_samples <- fast_rmvnorm(1e5, twostage_map, laplace_cut_res$asympvar)

#### Plot Ellipses for the two stage approach ####

# ## plot 95% credible region using stat_ellipse
# df_samples <- data.frame(theta1 = cut_samples[,1], theta2 = cut_samples[,2])
# ggplot(df_samples, aes(x = theta1, y = theta2)) + theme_bw()  +
#   # geom_point(alpha = 0.3) +
#   stat_ellipse(level = 0.95, color = 'blue', type = 'norm', linewidth = 2) +
#   stat_ellipse(level = 0.95, type = 'norm', mapping = aes(x = theta1, y = theta2),
#                data = data.frame(theta1 = cut_laplace_samples[,1], theta2 = cut_laplace_samples[,2]),
#                inherit.aes = FALSE,
#                color = 'red', linetype = 2, linewidth = 2) +
#   xlab(expression(theta[1])) + ylab(expression(theta[2])) +
#   ggtitle('95% credible regions: Cut Posterior (blue) vs Cut-Laplace (red)')
#
#### PBMI ####

pbmi_nsamples <- 5000
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("biased_dataset"))
clusterEvalQ(cl, {
  library(asympcuts)
})
pbmi_samples <- foreach (isample = 1:pbmi_nsamples, .combine = rbind) %dorng% {
  attach(biased_dataset)
  w_0 <- 1
  w1 <- rexp(n1, rate = 1)
  f_1 <- function(theta1) -sum(w1 * dnorm(x1, mean = theta1, sd = 1, log = TRUE)) - w_0 * dnorm(theta1, mean = 0, sd = sigma_1_prior, log = TRUE)
  optim_1 <- tryCatch(optim(par = 0, fn = f_1, method = 'BFGS'), error = function(e)  return(list(error = TRUE)))
  pbmi_theta1 <- optim_1$par
  # refresh weights
  w2 <- rexp(n2, rate = 1)
  f_2 <- function(theta2) -sum(w2 * dnorm(x2, mean = pbmi_theta1 + theta2, sd = 1, log = TRUE)) - w_0 * dnorm(theta2, mean = 0, sd = sigma_2_prior, log = TRUE)
  optim_2 <- tryCatch(optim(par = 0, fn = f_2, method = 'BFGS'), error = function(e)  return(list(error = TRUE)))
  pbmi_theta2 <- optim_2$par
  detach(biased_dataset)
  c(pbmi_theta1, pbmi_theta2)
}
#
# end cluster
stopCluster(cl)
#
## plot 95% credible region associated with all_samples using stat_ellipse
# df_samples <- data.frame(theta1 = cut_samples[,1], theta2 = cut_samples[,2])
# ggplot(df_samples, aes(x = theta1, y = theta2)) + theme_bw()  +
#   # geom_point(alpha = 0.3) +
#   stat_ellipse(level = 0.95, color = 'blue', type = 'norm', linewidth = 2) +
#   stat_ellipse(level = 0.95, type = 'norm', mapping = aes(x = theta1, y = theta2),
#                data = data.frame(theta1 = pbmi_samples[,1], theta2 = pbmi_samples[,2]),
#                inherit.aes = FALSE,
#                color = 'green', linetype = 2, linewidth = 2) +
#   stat_ellipse(level = 0.95, type = 'norm', mapping = aes(x = theta1, y = theta2),
#                data = data.frame(theta1 = cut_laplace_samples[,1], theta2 = cut_laplace_samples[,2]),
#                inherit.aes = FALSE,
#                color = 'red', linetype = 2, linewidth = 2) +
#   xlab(expression(theta[1])) + ylab(expression(theta[2])) +
#   ggtitle('95% credible regions: Cut Posterior (blue) vs PBMI (green) vs Cut-Laplace (red)')

laplace_cut_var <- laplace_cut_res$asympvar
save(cut_samples, twostage_map, laplace_cut_var, pbmi_samples, file = paste0("inst/biased_data/biased_data_v2_twostage_n1_", n1, ".RData"))
detach(biased_dataset)

