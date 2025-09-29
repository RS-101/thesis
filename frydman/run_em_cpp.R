# optional: quick check
mdl_ptr <- setup_frydman_cpp(res_w_ic$obs)
my_model <- model_data_summary(mdl_ptr)
setup <- model_data_to_list(mdl_ptr)

# Optionally provide initials; otherwise defaults inside C++ are used.
fit <- em_fit(mdl_ptr, z_init = rep(1/my_model$I_mark, my_model$I_mark), lambda_init = rep(0.1, length(my_model$N)),
              max_iter = 100, tol = 1e-8, verbose = TRUE)

str(fit, 1)
fit$z_i
fit$lambda_n

estimators <- calc_F_and_hazards(
  grid_points = seq(5, 15, by = 0.01),
  z_i = fit$z_i, lambda = fit$lambda_n, setup$Q, setup$Q_i_mark , setup$T_star, setup$E_star
)

plot_estimators_gg(estimators)

list(setup$Q, setup$Q_i_mark)




