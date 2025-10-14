source("frydman/helper_function_cpp.R")
source("frydman/helper_functions.R")
sourceCpp("frydman/functions_em.cpp")

get_npmle <- function(data, max_iter = 200, tol = 1e-8, verbose = FALSE) {
  
  mdl_ptr <- setup_frydman_cpp(data)
  my_model <- model_data_summary(mdl_ptr)
  setup <- model_data_to_list(mdl_ptr)
  
  # Optionally provide initials; otherwise defaults inside C++ are used.
  fit <- em_fit(mdl_ptr,
                z_init = rep(1/my_model$I_mark, my_model$I_mark),
                lambda_init = rep(0.1, length(my_model$N)),
                max_iter = max_iter, 
                tol = tol, 
                verbose = verbose)
  
  estimators <- calc_F_and_hazards(
    grid_points = seq(0, max(data$T_obs), length.out = 256),
    z_i = fit$z_i,
    lambda = fit$lambda_n, 
    setup$Q,
    setup$Q_i_mark, 
    setup$T_star, 
    setup$E_star,
    as_function = TRUE
  )
  
  list(
    estimators = estimators,
    settings = list(
      data = data,
      verbose = verbose,
      max_iter = max_iter,
      tol = tol
    ))
}





