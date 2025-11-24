library(asympcuts)
rm(list = ls())
#### Description of model ####
# x_{1,i} ~ Normal(\theta_1, \sigma_1^2), i = 1,...,n_1
# x_{2,i} ~ Normal(\theta_1 + \theta_2, \sigma_2^2), i = 1,...,n_2
# \theta_1 ~ Normal(0,1)
# \theta_2 ~ Normal(0,.1^2)
# n_1 = 20 and then 100
# n_2 = 10 * n_1 = 200 and then 1000


###
set.seed(1)
n1 <- 20
n2 <- 10 * n1
x1 <- rnorm(n1, 0, 1)
x2 <- rnorm(n2, 1, 0.5)
sigma_1_prior <- 1
sigma_2_prior <- 0.1
biased_dataset <- list(n1 = n1, n2 = n2, x1 = x1, x2 = x2, sigma_1_prior = sigma_1_prior, sigma_2_prior = sigma_2_prior)
save(biased_dataset, file = paste0("inst/biased_data/biased_data_n1_", n1, ".RData"))
rm(list = setdiff(ls(), c("biased_dataset")))
# rm(list = c("n1", "n2", "x1", "x2", "sigma_1_prior", "sigma_2_prior", "biased_dataset"))
## remove everything apart from biased_dataset

###
## Run standard Bayesian inference + Laplace approximation
## this creates inst/biased_data/biased_data_v2_jointmodel.RData
source("inst/biased_data/biased_data_v2_jointmodel.R")
rm(list = setdiff(ls(), c("biased_dataset")))
## Run two-stage Bayesian inference + Laplace approximation + PBMI
source("inst/biased_data/biased_data_v2_twostage.R")
## this creates inst/biased_data/biased_data_v2_twostage.RData
## Plots
source("inst/biased_data/biased_data_v2_plot.R")

set.seed(1)
n1 <- 100
n2 <- 10 * n1
x1 <- rnorm(n1, 0, 1)
x2 <- rnorm(n2, 1, 0.5)
sigma_1_prior <- 1
sigma_2_prior <- 0.1
biased_dataset <- list(n1 = n1, n2 = n2, x1 = x1, x2 = x2, sigma_1_prior = sigma_1_prior, sigma_2_prior = sigma_2_prior)
save(biased_dataset, file = paste0("inst/biased_data/biased_data_n1_", n1, ".RData"))
rm(list = setdiff(ls(), c("biased_dataset")))
## Run standard Bayesian inference + Laplace approximation
## this creates inst/biased_data/biased_data_v2_jointmodel.RData
source("inst/biased_data/biased_data_v2_jointmodel.R")
rm(list = setdiff(ls(), c("biased_dataset")))
## run two-stage inference (cut posterior, its Laplace approximation, and PBMI)
source("inst/biased_data/biased_data_v2_twostage.R")
## this creates inst/biased_data/biased_data_v2_twostage.RData
## Plots
source("inst/biased_data/biased_data_v2_plot.R")




