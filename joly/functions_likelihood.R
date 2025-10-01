library(splines2)
#### Create spline hazard ####
# needs testing
make_spline_funs <- function(param, knots, degree) {
  
  spline_hazard <- function(x, theta, knots, degree = 3, cumulative = F) {
    n_knots <- length(knots)
    if (length(theta) != degree + n_knots - 1) stop("theta wrong length")
    
    if (cumulative) {
      I <- iSpline(x,
                   degree = degree,
                   knots = knots[-c(1,n_knots)],
                   Boundary.knots = knots[c(1,n_knots)],
                   intercept = T,
                   warn.outside = F)
      return(as.vector(I %*% theta))
    } else {
      M <- mSpline(x,
                   degree = degree,
                   knots = knots[-c(1,n_knots)],
                   Boundary.knots = knots[c(1,n_knots)],
                   intercept = T,
                   warn.outside = F)
      return(as.vector(M %*% theta))
    }
  }
  
  
  
  return_obj <- list(
    hazard = function(x) spline_hazard(x, param, knots, degree, cumulative = FALSE),
    cumhaz = function(x) spline_hazard(x, param, knots, degree, cumulative = TRUE)
  )
  
  
  class(return_obj) <- c(class(return_obj), "single_spline_hazard")
  return_obj
}


# we assume order a12, a13 and a23
make_all_3_spline_funs <- function(theta, knots, degree) {
  # define all transitions
  a12_funs <- make_spline_funs(theta$a12, knots$a12, degree)
  a13_funs <- make_spline_funs(theta$a13, knots$a13, degree)
  a23_funs <- make_spline_funs(theta$a23, knots$a23, degree)
  
  # hazards
  a12 <- a12_funs$hazard
  a13 <- a13_funs$hazard
  a23 <- a23_funs$hazard
  
  # cumulative hazards
  A12 <- a12_funs$cumhaz
  A13 <- a13_funs$cumhaz
  A23 <- a23_funs$cumhaz
  
  return_obj <- list(a12 = a12,
                     a13 = a13, 
                     a23 = a23,
                     A12 = A12,
                     A13 = A13,
                     A23 = A23)
  
  class(return_obj) <- c(class(return_obj), "full_spline_hazard")
  return_obj
}

#### Likelihood functions for each case ####
case_1_likelihood <- function(V_0, V_healthy, T_obs, a12, a13, a23, A12, A13, A23) {
  if (all(V_healthy == T_obs)) {
    exp(-A12(V_healthy) - A13(V_healthy)) / exp(-A12(V_0) - A13(V_0))
  } else {
    p1 <- 1 / exp(-A12(V_0) - A13(V_0))
    p2 <- exp(-A12(T_obs) - A13(T_obs))
    
    factored_out <- exp(-A23(T_obs))
#    grid <- sort(unique(c(V_0, V_healthy, T_obs)))
    grid <- seq(min(V_healthy), max(T_obs), length.out = 250)
      
      
    integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))
    
    dx <- diff(grid)
    mid <- (integrand[-1] + integrand[-length(grid)]) / 2
    integral_at_grid <- c(0, cumsum(dx * mid))
    
    integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
    
    p3 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))
    
    p1 * (p2 + p3)
  }
}
case_2_likelihood <- function(V_0, V_healthy, T_obs, a12, a13, a23, A12, A13, A23) {
  p1 <- 1 / exp(-A12(V_0) - A13(V_0))
  p2 <- exp(-A12(T_obs) - A13(T_obs)) * a13(T_obs)
  
  factored_out <- exp(-A23(T_obs)) * a23(T_obs)
  #grid <- sort(unique(c(V_0, V_healthy, T_obs)))
  grid <- seq(min(V_healthy), max(T_obs), length.out = 250)
  integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))
  
  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))
  
  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
  
  p3 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))
  
  p1 * (p2 + p3)
}
case_3_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a12, a13, a23, A12, A13, A23) {
  p1 <- 1 / exp(-A12(V_0) - A13(V_0))
  
  factored_out <- exp(-A23(T_obs))
  #grid <- sort(unique(c(V_0, V_healthy, V_ill, T_obs)))
  grid <- seq(min(V_healthy), max(V_ill), length.out = 250)
  integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))
  
  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))
  
  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
  
  
  
  p2 <- factored_out * (integral_fun(V_ill) - integral_fun(V_healthy))
  
  p1 * p2
}
case_4_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a12, a13, a23, A12, A13, A23) {
  p1 <- 1 / exp(-A12(V_0) - A13(V_0))
  
  ####
  
  factored_out <- exp(-A23(T_obs)) * a23(T_obs)
  grid <- seq(min(V_healthy), max(V_ill), length.out = 250)
  integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))
  
  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))
  
  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
  
  
  
  p2 <- factored_out * (integral_fun(V_ill) - integral_fun(V_healthy))
  
  p1 * p2
}


#### Penalizing ####

# needs testing
pen_mat_m_splines <- function(input_knots, degree = 3) {
  n_knots <- length(input_knots)
  d <- diff(input_knots)
  g_ab <- mSpline(x = input_knots,
                  knots = input_knots[-c(1,n_knots)],
                  Boundary.knots = range(input_knots),
                  degree = degree,
                  derivs = 2,
                  intercept = TRUE
  )
  knots_mid <- input_knots[-length(input_knots)] + d / 2
  g_ab_mid <- mSpline(x = knots_mid,
                      knots = input_knots[-c(1,n_knots)],
                      Boundary.knots = range(input_knots),
                      degree = degree,
                      derivs = 2,
                      intercept = TRUE
  )
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}

#### Full log likelihood ####
# needs testing
# V_healthy = V_k or V_m
# V_ill = V_k+1


full_log_likehood <- function(V_0,
                              V_healthy, 
                              V_ill, 
                              T_obs, 
                              status, 
                              theta, 
                              degree,
                              knots) {
  if (!(length(V_0) == length(V_healthy) & length(V_healthy) == length(T_obs) & length(T_obs) == length(status))) stop("wrong lengths")

  # define all transitions
  a12_funs <- make_spline_funs(theta$a12, knots$a12, degree)
  a13_funs <- make_spline_funs(theta$a13, knots$a13, degree)
  a23_funs <- make_spline_funs(theta$a12, knots$a12, degree)

  # hazards
  a12 <- a12_funs$hazard
  a13 <- a13_funs$hazard
  a23 <- a23_funs$hazard

  # cumulative hazards
  A12 <- a12_funs$cumhaz
  A13 <- a13_funs$cumhaz
  A23 <- a23_funs$cumhaz
  
  loglik <- 0
  if (any(status == 1)) loglik <- loglik + sum(log(case_1_likelihood(V_0[as.integer(status) == 1],
                                                                     V_healthy[as.integer(status) == 1],
                                                                     T_obs[as.integer(status) == 1],
                                                                     a12, a13, a23, 
                                                                     A12, A13, A23)))
  
  if (any(as.integer(status) == 2)) loglik <- loglik + sum(log(case_2_likelihood(V_0[as.integer(status) == 2],
                                                                     V_healthy[as.integer(status) == 2],
                                                                     T_obs[as.integer(status) == 2],
                                                                     a12, a13, a23, 
                                                                     A12, A13, A23)))
  
  if (any(as.integer(status) == 3)) loglik <- loglik + sum(log(case_3_likelihood(V_0[as.integer(status) == 3],
                                                                     V_healthy[as.integer(status) == 3],
                                                                     V_ill[as.integer(status) == 3],
                                                                     T_obs[as.integer(status) == 3],
                                                                     a12, a13, a23, 
                                                                     A12, A13, A23)))
  
  if (any(as.integer(status) == 4)) loglik <- loglik + sum(log(case_4_likelihood(V_0[as.integer(status) == 4],
                                                                     V_healthy[as.integer(status) == 4],
                                                                     V_ill[as.integer(status) == 4],
                                                                     T_obs[as.integer(status) == 4],
                                                                     a12, a13, a23, 
                                                                     A12, A13, A23)))
  
  hazards = list(a12 = a12,
                 a13 = a13,
                 a12 = a23,
                 A12 = A12,
                 A13 = A13,
                 A12 = A23)
  
  class(hazards) <- c(class(hazards), "full_spline_hazard")

  return(list(loglik = loglik,
              theta = theta,
              knots = knots, 
              degree = degree,
              hazards = hazards))
}




#### Fit with penalty using base::optim ####
do_likelihood_optim <- function(sim_dat, n_knots, degree = 3, penalizer, knots = NULL) {
  n_par <- n_knots + degree - 1
  stopifnot(length(penalizer) == 3)
  
  # Use training-derived knots (important for CV fairness). Allow override.
  if (is.null(knots)) {
    knots_a12 <- seq(min(sim_dat$V_0),  max(sim_dat$T_obs), length.out = n_knots)
    knots_a13 <- knots_a12
    knots_a23 <- knots_a12
    knots <- list(a12 = knots_a12, a13 = knots_a13, a23 = knots_a23)
  } else {
    stopifnot(length(knots$a12) == n_knots, length(knots$a13) == n_knots, length(knots$a23) == n_knots)
  }
  
  # penalty matrix for each hazard (must exist in your env)
  P12 <- pen_mat_m_splines(knots$a12)
  P13 <- pen_mat_m_splines(knots$a13)
  P23 <- pen_mat_m_splines(knots$a23)
  
  # NOTE: We MINIMIZE:  -loglik + lambda * penalty_sum
  # Your original code subtracted the penalty and also used '-' between hazards,
  # which effectively rewards roughness. We fix both issues here.
  obj_fun <- function(x) {
    theta <- list(
      a12 = x[1:n_par],
      a13 = x[(n_par+1):(2*n_par)],
      a23 = x[(2*n_par+1):(3*n_par)]
    )
    ll <- full_log_likehood(
      V_0       = sim_dat$V_0,
      V_healthy = sim_dat$V_healthy,
      V_ill     = sim_dat$V_ill,
      T_obs     = sim_dat$T_obs,
      status    = sim_dat$status,
      theta     = theta,
      degree    = degree,
      knots     = knots
    )$loglik
    
    pen <- c(drop(t(theta$a12) %*% P12 %*% theta$a12),
             drop(t(theta$a13) %*% P13 %*% theta$a13),
             drop(t(theta$a23) %*% P23 %*% theta$a23))
    
    val <- -(ll - sum(penalizer * pen))
    if (!is.finite(val)) return(1e10)
    val
  }
  
  x0 <- rep(1, 3*n_par)
  # Warm check
  invisible(obj_fun(x0))
  
  res <- optim(
    par   = x0,
    fn    = obj_fun,
    method= "L-BFGS-B",
    lower = rep(0.0001, 3*n_par),
    upper = rep(13,     3*n_par)
  )
  
  theta_hat <- list(
    a12 = res$par[1:n_par],
    a13 = res$par[(n_par+1):(2*n_par)],
    a23 = res$par[(2*n_par+1):(3*n_par)]
  )
  
  pen_term_at_theta_hat = sum(penalizer * c(drop(t(theta_hat$a12) %*% P12 %*% theta_hat$a12),
                                            drop(t(theta_hat$a13) %*% P13 %*% theta_hat$a13),
                                            drop(t(theta_hat$a23) %*% P23 %*% theta_hat$a23)))

  est_hazards <- make_all_3_spline_funs(theta_hat, knots, degree)
  
  list(
    loglik = -res$value + pen_term_at_theta_hat,
    pen_loglik = -res$value,
    hazards = est_hazards,
    theta        = theta_hat,
    knots        = knots,
    degree       = degree
  )
}
