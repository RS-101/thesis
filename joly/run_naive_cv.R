source("datageneration/functions_simulate_data.R")
debugSource("joly/function_perform_cv.R")

library(tidyverse)

simulated_data <- simulate_idm_weibull(1000)

# Choose your grids (adjust as needed)
# n_knots must be >= 3; include a range to let CV choose smoothness via basis size + lambda penalizer grid on log-scale is typical
grid_nk <- c(5)
grid_lam <- 10^seq(3,4, length.out = 2)



set.seed(21321)
cv_res <- cv_penalized_splines(simulated_data$obs, 
                               n_knots_grid = grid_nk,
                               penalizer_a12 = grid_lam,
                               penalizer_a13 = grid_lam,
                               penalizer_a23 = grid_lam,
                               K = 5)

cv_res$cv_table
cv_res$best

res_full <- cv_res$final_fit


