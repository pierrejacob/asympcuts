library(asympcuts)

rm(list = ls())

#### Description of model ####
## z_i ~ Normal(theta_1, 1), i=1,..,n
## y_i ~ Normal(theta_1 + theta_2 x_i, 1), i=1,..,n
## Priors: Normal(0, 10^2) on theta1 and theta2
## DGP: x_1,...,x_n ~ Normal(1,1), and (x_i,y_i) ~ Normal((0, 1 + x_i), Sigma) with
## Sigma = ( 1    rho )
##         ( rho  sigma^2)
## and we consider two scenarios:
## SCENARIO 1 : rho = 0.0, sigma = 1,
## SCENARIO 2 : rho = 0.5, sigma = 2.

## functions implementing the log-likelihoods and log-priors
ll1 <- function(theta1, z){
  sum(dnorm(z, mean = theta1, sd = 1, log = TRUE))
}
ll2 <- function(theta1, theta2, x, y){
  sum(dnorm(y, mean = theta1 + theta2 * x, sd = 1, log = TRUE))
}
lp1 <- function(theta1){
  dnorm(theta1, mean = 0, sd = 10, log = TRUE)
}
lp2 <- function(theta2){
  dnorm(theta2, mean = 0, sd = 10, log = TRUE)
}


## set seed and generate data under two scenarios
## this creates inst/toy_example/toy_example_data.RData
source("inst/toy_example/toy_example_v2_generate_data.R")
## remove objects from environment apart from ll1, ll2, lp1, lp2, and SCENARIO
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2")))
##
SCENARIO <- 1
## run joint model inference (standard posterior and its Laplace approximation)
## this creates inst/toy_example/toy_example_v2_jointmodel_scenario1.RData
source("inst/toy_example/toy_example_v2_jointmodel.R")
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2", "SCENARIO")))
## run two-stage inference (cut posterior, its Laplace approximation, and PBMI)
## this creates inst/toy_example/toy_example_v2_twostage_scenario1.RData
source("inst/toy_example/toy_example_v2_twostage.R")
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2", "SCENARIO")))
## run coverage experiments
## this creates inst/toy_example/toy_example_v2_twostage_coverage_scenario1.RData
source("inst/toy_example/toy_example_v2_twostage_coverage.R")
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2", "SCENARIO")))
## generate plots and tables
source("inst/toy_example/toy_example_v2_plot.R")

SCENARIO <- 2
## run joint model inference (standard posterior and its Laplace approximation)
## this creates inst/toy_example/toy_example_v2_jointmodel_scenario2.RData
source("inst/toy_example/toy_example_v2_jointmodel.R")
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2", "SCENARIO")))
## run two-stage inference (cut posterior, its Laplace approximation, and PBMI)
## this creates inst/toy_example/toy_example_v2_twostage_scenario2.RData
source("inst/toy_example/toy_example_v2_twostage.R")
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2", "SCENARIO")))
## run coverage experiments
## this creates inst/toy_example/toy_example_v2_twostage_coverage_scenario2.RData
source("inst/toy_example/toy_example_v2_twostage_coverage.R")
rm(list = setdiff(ls(), c("ll1", "ll2", "lp1", "lp2", "SCENARIO")))
## generate plots and tables
source("inst/toy_example/toy_example_v2_plot.R")




