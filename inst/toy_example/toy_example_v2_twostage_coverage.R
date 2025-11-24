## number of credible regions generated to compute coverage
ncrs <- 1000


set.seed(1)

# sample size
n1 <- n2 <- n <- 1000
load("inst/toy_example/toy_example_data.RData")
if (SCENARIO == 1){
  ##### rho = 0 and sigma^2 = 1
  rho <- 0.0
  sigma <- 1
}
if(SCENARIO == 2){
  ##### rho = 0.5 and sigma^2 = 2
  rho <- 0.5
  sigma <- sqrt(2)
}

##### Limiting values of the 2SMAP #####
theta1star <- meanz <- 0
# recall x ~ Normal(1,1)
meanx <- 1
varx <- 1
meanx2 <- varx + meanx^2
# recall y = 1 + x + Normal(0,1)
# so xy = (1 + x + Normal(0,1))*x = x + x^2 + x*Normal(0,1)
meany <- 2
meanxy <- 3
theta2star <- (meanxy - meanx * meanz) / meanx2
thetastar <- c(theta1star, theta2star)


## generate data
gendata <- function(n, rho, sigma){
  z_y_mat <- mvtnorm::rmvnorm(n, c(0,0), matrix(c(1, rho, rho, sigma^2), nrow =2))
  z <- meanz + z_y_mat[1:n1,1]
  x <- meanx + sqrt(varx) * rnorm(n, 0, 1)
  y <-  x + z_y_mat[1:n,2] + 1
  return(list(z = z, x = x, y = y))
}

#### Cut Laplace approximation ####

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, "SCENARIO")
clusterEvalQ(cl, {
  library(asympcuts)
})

credible_regions_laplace <- foreach(icr = 1:ncrs, .combine = rbind) %dorng% {
  data_test <- gendata(n, rho, sigma)
  x <- data_test$x
  y <- data_test$y
  z <- data_test$z
  laplace_cut_res <- laplace_cut(logposterior1 = function(theta1) ll1(theta1, z) + lp1(theta1),
                                 logposterior2 = function(theta1, theta2) ll2(theta1, theta2, x, y) + lp2(theta2),
                                 init_theta1 = 0,
                                 init_theta2 = 0)
  twostage_map <- laplace_cut_res$thetahat
  cut_cov_laplace <- laplace_cut_res$asympvar
  ## exact CR for theta1star
  cr1 <- c(qnorm(0.025, mean = twostage_map[1], sd = sqrt(cut_cov_laplace[1,1])), qnorm(0.975, mean = twostage_map[1], sd = sqrt(cut_cov_laplace[1,1])))
  ## exact CR for theta2star
  cr2 <- c(qnorm(0.025, mean = twostage_map[2], sd = sqrt(cut_cov_laplace[2,2])), qnorm(0.975, mean = twostage_map[2], sd = sqrt(cut_cov_laplace[2,2])))
  ## exact CR for theta1star + theta2star
  var_sum_theta <- cut_cov_laplace[1,1] + cut_cov_laplace[2,2] + 2 * cut_cov_laplace[1,2]
  cr1p2 <- c(qnorm(0.025, mean = sum(twostage_map), sd = sqrt(var_sum_theta)), qnorm(0.975, mean = sum(twostage_map), sd = sqrt(var_sum_theta)))
  df_ <- data.frame(rbind(cr1, cr2, cr1p2))
  colnames(df_) <- c("lower", "upper")
  df_$parameter <- c("theta1", "theta2", "theta1p2")
  df_$icr <- icr
  df_
}
stopCluster(cl)
row.names(credible_regions_laplace) <- NULL
head(credible_regions_laplace)


## Compute coverage

## For theta1:
credible_regions_laplace %>% filter(parameter == "theta1") %>%
  mutate(coverage = ifelse((lower <= theta1star) & (upper >= theta1star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage))
## For theta2:
credible_regions_laplace %>% filter(parameter == "theta2") %>%
  mutate(coverage = ifelse((lower <= theta2star) & (upper >= theta2star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage))
## For theta1 + theta2:
credible_regions_laplace %>% filter(parameter == "theta1p2") %>%
  mutate(coverage = ifelse((lower <= (theta1star + theta2star)) & (upper >= (theta1star + theta2star)), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage))

#### PBMI ####

pbmi_nsamples <- 1000

cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, "SCENARIO")
clusterEvalQ(cl, {
  library(asympcuts)
  library(numDeriv)
})

credible_regions_pbmi <- foreach(icr = 1:ncrs, .combine = rbind) %dorng% {
  data_test <- gendata(n, rho, sigma)
  x <- data_test$x
  y <- data_test$y
  z <- data_test$z

  pbmi_samples <- foreach (isample = 1:pbmi_nsamples, .combine = rbind) %do% {
    w_0 <- 1
    w1 <- rexp(n, rate = 1)
    f_1 <- function(theta1) -sum(w1 * dnorm(z, mean = theta1, sd = 1, log = TRUE)) - w_0 * dnorm(theta1, mean = 0, sd = 10, log = TRUE)
    optim_1 <- tryCatch(optim(par = 0, fn = f_1, method = 'BFGS'), error = function(e)  return(list(error = TRUE)))
    pbmi_theta1 <- optim_1$par
    # refresh weights
    # w1 <- rexp(n, rate = 1)
    f_2 <- function(theta2) -sum(w1 * dnorm(y, mean = pbmi_theta1 + theta2 * x, sd = 1, log = TRUE)) - w_0 * dnorm(theta2, mean = 0, sd = 10, log = TRUE)
    optim_2 <- tryCatch(optim(par = 0, fn = f_2, method = 'BFGS'), error = function(e)  return(list(error = TRUE)))
    pbmi_theta2 <- optim_2$par
    c(pbmi_theta1, pbmi_theta2)
  }
  cr1 <- quantile(pbmi_samples[,1], probs = c(0.025, 0.975))
  cr2 <- quantile(pbmi_samples[,2], probs = c(0.025, 0.975))
  cr1p2 <- quantile(pbmi_samples[,1] + pbmi_samples[,2], probs = c(0.025, 0.975))
  df_ <- data.frame(rbind(cr1, cr2, cr1p2))
  colnames(df_) <- c("lower", "upper")
  df_$parameter <- c("theta1", "theta2", "theta1p2")
  df_$icr <- icr
  df_
}
stopCluster(cl)
row.names(credible_regions_pbmi) <- NULL
head(credible_regions_pbmi)

## Compute coverage
## For theta1:
credible_regions_pbmi %>% filter(parameter == "theta1") %>%
  mutate(coverage = ifelse((lower <= theta1star) & (upper >= theta1star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage))
## For theta2:
credible_regions_pbmi %>% filter(parameter == "theta2") %>%
  mutate(coverage = ifelse((lower <= theta2star) & (upper >= theta2star), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage))
## For theta1 + theta2:
credible_regions_pbmi %>% filter(parameter == "theta1p2") %>%
  mutate(coverage = ifelse((lower <= (theta1star + theta2star)) & (upper >= (theta1star + theta2star)), 1, 0)) %>%
  summarise(coverage_prob = mean(coverage))


save(theta1star, theta2star, thetastar, SCENARIO, n,
     credible_regions_laplace,
     credible_regions_pbmi,
     pbmi_nsamples,
     file = paste0("inst/toy_example/toy_example_v2_twostage_coverage_scenario",
                   SCENARIO, ".RData"))

