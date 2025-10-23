source("frydman/helper_functions.R")
sourceCpp("optimizing_function_em.cpp")

get_npmle <- function(data, data_list = NULL, max_iter = 200, tol = 1e-8, verbose = FALSE) {
  
  if(is.null(data_list )){
    data_list <- setup_data_to_list_format(data)
  }
  mdl_ptr <- make_model_data(data_list)
  
  z_init <- runif(data_list$I_mark)
  
  z_init <- z_init/sum(z_init)
  
  lambda_init <- runif(data_list$N,min = 0.1, max = 0.5)
  
  # Optionally provide initials; otherwise defaults inside C++ are used.
  fit <- em_fit(mdl_ptr,
                z_init = z_init,
                lambda_init = lambda_init,
                max_iter = max_iter, 
                tol = tol, 
                verbose = verbose)
  
  estimators <- calc_F_and_hazards(
    grid_points = seq(0, max(data$T_obs), length.out = 256),
    z_i = fit$z_i,
    lambda = fit$lambda_n, 
    data_list$Q,
    data_list$Q_i_mark, 
    data_list$t_star_n, 
    data_list$E_star
  )
  
  list(
    estimators = estimators$as_functions,
    estimators_point = estimators$as_points,
    raw_em_res = list(
      z_i = fit$z_i,
      lambda = fit$lambda_n
    ),
    settings = list(
      data = data,
      verbose = verbose,
      max_iter = max_iter,
      tol = tol
    ))
}
