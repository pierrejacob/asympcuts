## Compute Laplace approximation to the cut distribution
# ## posterior mode
# theta1hat <- (posterior_phi_alpha - 1) / (posterior_phi_alpha + posterior_phi_beta - 2)
# ## Boundary problem for 10-th parameter

# ## posterior mean
theta1hat <- posterior_phi_alpha / (posterior_phi_alpha + posterior_phi_beta)
#
init_theta2 <- rnorm(2)

optim_stage2 <- tryCatch(optim(par = init_theta2,
                               fn = function(theta2) -logposterior2(theta1hat, theta2),
                               method = 'BFGS'),
                         error = function(e)  return(list(error = TRUE)))
theta2hat <- optim_stage2$par # second MAP
theta2hat

d1 <- length(theta1hat)
d2 <- length(theta2hat)

# load('inst/epidemio/epidemio_cut_mcmc_samples.RData')
# theta2hat <- colMeans(theta2s)

thetahat <- c(theta1hat, theta2hat)
## Computes hessians at the two stage MAP estimate
hessian1 <- hessian(func = logposterior1, thetahat[1:d1])
# rootSolve::hessian(logposterior1, thetahat[1:d1], centered = FALSE)
## Compute Hessian of L_2 at 2SMAP
hessian2 <- hessian(func = function(theta) logposterior2(theta[1:d1], theta[(d1+1):(d1+d2)]), thetahat)
  ## Construct asymptotic precision matrix of the Laplace approximation
asympprec <- matrix(0, nrow = d1 + d2, ncol = d1 + d2)
topleft <- hessian1 + hessian2[1:d1,(d1+1):(d1+d2),drop=F] %*% solve(hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2),drop=F], hessian2[(d1+1):(d1+d2),1:d1,drop=F])
asympprec[1:d1, 1:d1] <- topleft
asympprec[1:d1,(d1+1):(d1+d2)] <- hessian2[1:d1,(d1+1):(d1+d2)]
asympprec[(d1+1):(d1+d2),1:d1] <- hessian2[(d1+1):(d1+d2),1:d1]
asympprec[(d1+1):(d1+d2),(d1+1):(d1+d2)] <- hessian2[(d1+1):(d1+d2),(d1+1):(d1+d2)]
## Laplace approximation of cut posterior is given by N(thetahat, asympvar)
asympvar = solve(-asympprec)

save(thetahat, asympvar, file = "inst/epidemio/epidemio_cutlaplace.RData")

laplace_samples <- fast_rmvnorm(n = 1e5, thetahat[14:15], asympvar[(d1+1):(d1+d2),(d1+1):(d1+d2)])
# laplace_samples <- fast_rmvnorm(n = 1e5, colMeans(theta2s), asympvar[(d1+1):(d1+d2),(d1+1):(d1+d2)])
laplace_samples <- data.frame(theta1 = laplace_samples[,1], theta2 = laplace_samples[,2])

ggplot(laplace_samples, aes(x = theta1, y = theta2)) + theme_bw()  +
  stat_ellipse(level = 0.95, color = 'blue', type = 'norm', linewidth = 2) +
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  ggtitle('95% credible regions: Cut-Laplace Posterior (blue)')
#
