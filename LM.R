source("logger_function.R")

library(numDeriv)




# Gauss Newton algorithm
naive_GN <- function(f, x, y, p_0, lambda, max_rep) {
  
  logger <- create_logger()
  
  # f should be function(x, p)
  f_at_x = function(p) f(x, p)
  
  func_residual <- function(p) y-f_at_x(p)
  dim_p <- length(p_0)
  dim_x <- length(x)
  parameters <- vector("list", max_rep)
  parameters[[1]] = p_0
  
  
  for (i in seq_along(parameters)[-1]) {
    p <- parameters[[i-1]]
    
    calc_residual <- func_residual(p)
    J_at_p <- numDeriv::jacobian(f_at_x, p)

    N = crossprod(J_at_p) + diag(rep(lambda, dim_p))
    RHS <- crossprod(J_at_p, calc_residual)
    
    delta_p = solve(N, RHS)
    print(delta_p)
    parameters[[i]] = delta_p + p 
    
    
    logger$log(
      iter = i,
      loss = crossprod(calc_residual),
      extra = list(time = Sys.time()))
  }
  
  list(res = parameters, log = logger$get_log())
}


# Marquardt's algorithm
naive_LM <- function(f, x, y, p_0, lambda, max_rep) {
  logger <- create_logger()
  # f should be function(x, p)
  f_at_x = function(p) f(x, p)
  
  func_residual <- function(p) y-f_at_x(p)
  dim_p <- length(p_0)
  dim_x <- length(x)
  parameters <- vector("list", max_rep)
  parameters[[1]] = p_0
  
  
  for (i in seq_along(parameters)[-1]) {
    p <- parameters[[i-1]]
    
    calc_residual <- func_residual(p)
    size_calc_residual <- crossprod(calc_residual)
    J_at_p <- numDeriv::jacobian(f_at_x, p)

    N = crossprod(J_at_p)
    RHS <- crossprod(J_at_p, calc_residual)
    
    delta_p = solve(N + diag(rep(lambda, dim_p)), RHS)
    
    p_new = p + delta_p
    
    calc_residual_new = func_residual(p_new)
    size_calc_residual_new = crossprod(calc_residual_new)
    
    lambda_vec <- list()
    loss_vec <- list()
    j <- 1
    lambda_vec[[j]] = lambda
    loss_vec[[j]] = size_calc_residual
    
    while(size_calc_residual < size_calc_residual_new) {
      lambda <- lambda*1.5
      
      delta_p = solve(N + diag(rep(lambda, dim_p)), RHS)
      
      p_new = p + delta_p
      
      calc_residual_new = func_residual(p_new)
      size_calc_residual_new = crossprod(calc_residual_new)
      
      j = j + 1
      
      lambda_vec[[j]] = lambda
      loss_vec[[j]] = size_calc_residual_new
    }
    
    lambda <- lambda * 0.75
    parameters[[i]] = delta_p + p 
    
    logger$log(
      iter = i,
      loss = crossprod(calc_residual),
      extra = list(damping_loops = j,
                   lambda_vec = paste(unlist(lambda_vec), collapse = ", "),
                   loss_vec = paste(unlist(loss_vec), collapse = ", ")))
      }
  
  list(res = parameters, log = logger$get_log(), raw_log = logger$record)
}




