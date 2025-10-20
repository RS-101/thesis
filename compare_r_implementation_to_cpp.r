library(Rcpp)
sourceCpp("optimizing_function_em.cpp")

source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")
source("frydman/functions_em.R")
source("frydman/helper_functions.R")


sim_data <- simulate_idm_constant_hazards(n = 250, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016)

data_list <- setup_data_to_list_format(sim_data$obs, T)


set.seed(1)
z_init <- runif(data_list$cpp_data$I_mark)
z_init <- z_init/sum(z_init)

lambda_init <- rep(0.1,data_list$cpp_data$N)

max_iter = 1

res_cpp <- em_fit(make_model_data(data_list$cpp_data),
              z_init = z_init,
              lambda_init = lambda_init,
              max_iter = max_iter, 
              tol = 0, 
              verbose = T)

res_em <- do.call(
  em_estimate_raw,
  c(list(
    verbose = T,
    max_iter = max_iter,
    tol = 0,
    z_init = z_init,
    lambda_init = lambda_init
  ),data_list$r_data)
)

# change digits to see more detailed differences
options(digits = 16)
sum(res_em$z - res_cpp$z_i)

sum(res_em$z)
res_cpp$z_i


sum(abs(as.numeric(res_em$alpha_ij) - as.numeric(res_cpp$alpha_ij)))
sum(abs(as.numeric(res_em$beta_im) - as.numeric(res_cpp$beta_im)))
sum(abs(as.numeric(res_em$mu_mi) - as.numeric(res_cpp$mu_mi)))
sum(abs(as.numeric(res_em$eta_ui) - as.numeric(res_cpp$eta_ui)))
sum(abs(as.numeric(res_em$gamma_ci) - as.numeric(res_cpp$gamma_ci)))
