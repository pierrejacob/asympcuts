#'@rdname laplace_cut
#'@title Laplace approximation of the cut posterior
#'@description Centered at the MAP and more or less using the Hessian of the log-density of the cut posterior
#' @param logposterior1 : a function of theta1
#' @param logposterior2 : a function of theta1 and theta2
#' @param init_theta1 : initial point for theta1
#' @param init_theta2 : initial point for theta2
#' @return a list containing:
#'\itemize{
#' \item thetahat : the MAP estimate
#' \item asympvar : the asymptotic variance matrix of the Laplace approximation
#'}
#'@export
laplace_cut <- function(logposterior1, logposterior2, init_theta1, init_theta2){
  d1 <- length(init_theta1)
  d2 <- length(init_theta2)
  ##### Two-stage MAP estimator
  optim_stage1 <- tryCatch(optim(par = init_theta1,
                                 fn = function(theta1) -logposterior1(theta1),
                                 method = 'BFGS'),
                           error = function(e)  return(list(error = TRUE)))
  theta1hat <- optim_stage1$par # first MAP
  optim_stage2 <- tryCatch(optim(par = init_theta2,
                                 fn = function(theta2) -logposterior2(theta1hat, theta2),
                                 method = 'BFGS'),
                           error = function(e)  return(list(error = TRUE)))
  theta2hat <- optim_stage2$par # second MAP
  thetahat <- c(theta1hat, theta2hat)
  ## Computes hessians at the two stage MAP estimate
  hessian1 <- hessian(func = logposterior1, thetahat[1:d1])
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
  return(list(thetahat = thetahat, asympvar = solve(-asympprec)))
}
