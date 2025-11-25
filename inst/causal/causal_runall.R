library(asympcuts)
rm(list = ls())
set.seed(1)

## Run Cut-Posterior, Cut-Laplace and PBMI
## creates an .rds file and an RData file
source("inst/causal/causal_twostage.R")
## create all plots
source("inst/causal/causal_plot.R")

## optionally, create more plots on a synthetic experiment, in a different model
## where the propensity score are used as covariate in the second module,
## and where the Laplace approximation is accurate

# source("inst/causal/causal_synthetic_psascovariate.R")
