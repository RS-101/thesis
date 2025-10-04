library(Rcpp)
source("datageneration/functions_simulate_data.R")
sourceCpp("frydman/functions_em.cpp")
source("frydman/setup_cpp_from_data.R")
source("frydman/functions_extract_estimators.R")
source("frydman/functions_plot.R")
source("frydman/setup_for_R_from_data.R")

set.seed(100)
sim_dat <- simulate_idm_weibull(1000,
                                shape12 = 1.1, scale12 = 1/0.0008,
                                shape13 = 1.8, scale13 = 1/0.0002,
                                shape23 = 1.3, scale23 = 1/0.0016)

sim_dat$true_data_generation$plot

# optional: quick check
mdl_ptr <- setup_frydman_cpp(sim_dat$obs)
my_model <- model_data_summary(mdl_ptr)
setup <- model_data_to_list(mdl_ptr)

# Optionally provide initials; otherwise defaults inside C++ are used.
fit <- em_fit(mdl_ptr, z_init = rep(1/my_model$I_mark, my_model$I_mark), lambda_init = rep(0.1, length(my_model$N)),
              max_iter = 1000, tol = 0.00001, verbose = TRUE)

str(fit, 1)
fit$z_i
fit$lambda_n

estimators <- calc_F_and_hazards(
  grid_points = seq(0, 5000, length.out = 256),
  z_i = fit$z_i, lambda = fit$lambda_n, setup$Q, setup$Q_i_mark , setup$T_star, setup$E_star
)





plot_estimators_gg(estimators)

list(setup$Q, setup$Q_i_mark)


plot_cumhazards_est_vs_true(estimators, sim_dat$true_data_generation$hazards)

