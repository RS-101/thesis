library(qs)
library(Rcpp)

source("frydman/wrapper.R")
source("datageneration/functions_simulate_data.R")
source("experiments/plot_function.R")

iter <- 100

bias12 <- numeric(iter)
bias13  <- numeric(iter)
bias23 <- numeric(iter)


for (i in 1:iter) {
  sim_data <- simulate_idm_constant_hazards(n = 250, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016)
  A12_true <- integrate(sim_data$true_data_generation$hazards$a12, lower = 0, upper = 365)
  A13_true <- integrate(sim_data$true_data_generation$hazards$a13, lower = 0, upper = 365)
  A23_true <- integrate(sim_data$true_data_generation$hazards$a23, lower = 0, upper = 365)
  
  print(i)
  
  A12_est <- numeric(100)
  A13_est <- numeric(100)
  A23_est <- numeric(100)
  
  for(n in 1:100) {
    est_npmle <- get_npmle(sim_data$obs, 
                           max_iter = 100, tol = 0.03, verbose = FALSE)
    A12_est[n] <- est_npmle$estimators$Lambda12(365)
    A13_est[n] <- est_npmle$estimators$Lambda13(365)
    A23_est[n] <- est_npmle$estimators$Lambda23(365)
  }
  
  bias12[i] <- A12_true$value - mean(A12_est, na.rm = T)
  bias13[i] <- A13_true$value - mean(A13_est, na.rm = T)
  bias23[i] <- A23_true$value - mean(A23_est, na.rm = T)
}

mean(bias12)
mean(bias13)
mean(bias23)


