library(asympcuts)
library(Rcpp)
rm(list = ls())

# Plummer's example on HPV
# nhpv considered as Y
nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
y1 <- matrix(nhpv, nrow = 1)
# Npart is put in the parameters
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
           143, 229, 696, 93)
J <- 13
# posterior is beta in each study, with parameters
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv
## first log-posterior
logposterior1 <- function(theta1){
  sum(dbeta(theta1, shape1 = posterior_phi_alpha, shape2 = posterior_phi_beta, log = TRUE))
}

###
# For module 2, ncases considered data
ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
y2 <- matrix(ncases, nrow = 1)
# Npop considered parameters
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
          26751, 75815, 150302, 354993, 3683043, 507218)
Npop_normalized <- log(10**(-3) * Npop)
# Find parameters given Y2
hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))
dprior2 <- function(theta2, hyper2){
  return(sum(dnorm(theta2, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)))
}

sourceCpp('inst/epidemio/epidemiology.cpp')

# cppFunction("double plummer_module2_loglikelihood_(NumericVector theta1, NumericVector theta2,
#             NumericVector ncases, NumericVector Npop_normalized){
#     double eval = 0;
#     double mu, logmu;
#     for (int j = 0; j < 13; j++){
#       logmu = theta2(0) + theta1(j) * theta2(1) + Npop_normalized(j);
#       mu = exp(logmu);
#       eval += ncases(j) * logmu - mu;
#     }
#       return eval;
#     }
#       ")

logposterior2 <- function(theta1, theta2) epidemiology_loglikC(theta1, theta2, ncases, Npop_normalized) + dprior2(theta2, hyper2)


# Run one-stage MCMC to estimate the joint model's posterior distribution
# Generate epidemiology_joint_mcmc_samples.RData
source("inst/epidemio/epidemio_joint.R")
# Run two-stage MCMC to estimate the cut distribution
# Generate epidemiology_cut_mcmc_samples.RData
source("inst/epidemio/epidemio_twostage.R")
# Run PBMI
source("inst/epidemio/epidemio_pbmi.R")
# Run Cut-Laplace approximation
source("inst/epidemio/epidemio_cutlaplace.R")
# Generate plots
source("inst/epidemio/epidemio_plot.R")
