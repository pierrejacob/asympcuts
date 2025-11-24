library(asympcuts)
rm(list = ls())
set.seed(1)

## Run Cut-Posterior, Cut-Laplace and PBMI
## creates an .rds file and an RData file
source("inst/causal/causal_twostage.R")
## create all plots
source("inst/causal/causal_plot.R")

## optionally, create more plots on a synthetic experiment with more data, to see if Laplace approximation becomes more accurate
## takes quite a while as n is set to 25,000

# source("inst/causal/causal_synthetic_largen.R")
