library(asympcuts)
rm(list = ls())
set.seed(1)

## Run Cut-Posterior, Cut-Laplace and PBMI
## creates an .rds file and an RData file
source("inst/causal/causal_lalonde/causal_twostage.R")
## create all plots
source("inst/causal/causal_lalonde/causal_plot.R")
