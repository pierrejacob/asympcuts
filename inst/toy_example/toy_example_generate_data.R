
##### generate z and y #####
set.seed(1)
# sample size
n1 <- n2 <- n <- 1000
##### rho = 0 and sigma^2 = 1 #####
rho <- 0.0
sigma <- 1
## generate data
z_y_mat <- mvtnorm::rmvnorm(n, c(0,0), matrix(c(1, rho, rho, sigma^2), nrow =2))
meanz = 0
meanx = 1
varx = 1
z <- meanz + z_y_mat[1:n1,1]
x <- meanx + sqrt(varx) * rnorm(n, 0, 1)
y <-  x + z_y_mat[1:n,2] + 1

data_scenario1 <- list(n = n, rho = rho, sigma = sigma, x = x, y = y, z = z)
##### rho = 0.5 and sigma^2 = 2 #####
rho <- 0.5
sigma <- sqrt(2)
## generate data
z_y_mat <- mvtnorm::rmvnorm(n, c(0,0), matrix(c(1, rho, rho, sigma^2), nrow =2))
z <- meanz + z_y_mat[1:n1,1]
x <- meanx + sqrt(varx) * rnorm(n, 0, 1)
y <-  x + z_y_mat[1:n,2] + 1

data_scenario2 <- list(n = n, rho = rho, sigma = sigma, x = x, y = y, z = z)
save(data_scenario1, data_scenario2, file = "inst/toy_example/toy_example_data.RData")

