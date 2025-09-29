library(splines2)

case_1_likelihood <- function(V_0, V_healthy, T_obs, a01, a02, a12, A01, A02, A12) {
  if (all(V_healthy == T_obs)) {
    exp(-A01(V_healthy) - A02(V_healthy)) / exp(-A01(V_0) - A02(V_0))
  } else {
    p1 <- 1 / exp(-A01(V_0) - A02(V_0))
    p2 <- exp(-A01(T_obs) - A02(T_obs))
    
    factored_out <- exp(-A12(T_obs))
#    grid <- sort(unique(c(V_0, V_healthy, T_obs)))
    grid <- seq(min(V_healthy), max(T_obs), length.out = 250)
      
      
    integrand <- exp(-A01(grid) - A02(grid)) * a01(grid) * exp(A12(grid))
    
    dx <- diff(grid)
    mid <- (integrand[-1] + integrand[-length(grid)]) / 2
    integral_at_grid <- c(0, cumsum(dx * mid))
    
    integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
    
    p3 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))
    
    p1 * (p2 + p3)
  }
}
case_2_likelihood <- function(V_0, V_healthy, T_obs, a01, a02, a12, A01, A02, A12) {
  p1 <- 1 / exp(-A01(V_0) - A02(V_0))
  p2 <- exp(-A01(T_obs) - A02(T_obs)) * a02(T_obs)
  
  factored_out <- exp(-A12(T_obs)) * a12(T_obs)
  #grid <- sort(unique(c(V_0, V_healthy, T_obs)))
  grid <- seq(min(V_healthy), max(T_obs), length.out = 250)
  integrand <- exp(-A01(grid) - A02(grid)) * a01(grid) * exp(A12(grid))
  
  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))
  
  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
  
  p3 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))
  
  p1 * (p2 + p3)
}
case_3_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a01, a02, a12, A01, A02, A12) {
  p1 <- 1 / exp(-A01(V_0) - A02(V_0))
  
  factored_out <- exp(-A12(T_obs))
  #grid <- sort(unique(c(V_0, V_healthy, V_ill, T_obs)))
  grid <- seq(min(V_healthy), max(V_ill), length.out = 250)
  integrand <- exp(-A01(grid) - A02(grid)) * a01(grid) * exp(A12(grid))
  
  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))
  
  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
  
  
  
  p2 <- factored_out * (integral_fun(V_ill) - integral_fun(V_healthy))
  
  p1 * p2
}
case_4_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a01, a02, a12, A01, A02, A12) {
  p1 <- 1 / exp(-A01(V_0) - A02(V_0))
  
  ####
  
  factored_out <- exp(-A12(T_obs)) * a12(T_obs)
  grid <- seq(min(V_healthy), max(V_ill), length.out = 250)
  integrand <- exp(-A01(grid) - A02(grid)) * a01(grid) * exp(A12(grid))
  
  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))
  
  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")
  
  
  
  p2 <- factored_out * (integral_fun(V_ill) - integral_fun(V_healthy))
  
  p1 * p2
}

# needs testing
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

  make_spline_funs <- function(param, knots) {
    list(
      hazard = function(x) spline_hazard(x, param, knots, degree, cumulative = FALSE),
      cumhaz = function(x) spline_hazard(x, param, knots, degree, cumulative = TRUE)
    )
  }

  # define all transitions
  a01_funs <- make_spline_funs(theta$a01, knots$a01)
  a02_funs <- make_spline_funs(theta$a02, knots$a02)
  a12_funs <- make_spline_funs(theta$a12, knots$a12)

  # hazards
  a01 <- a01_funs$hazard
  a02 <- a02_funs$hazard
  a12 <- a12_funs$hazard

  # cumulative hazards
  A01 <- a01_funs$cumhaz
  A02 <- a02_funs$cumhaz
  A12 <- a12_funs$cumhaz
  
  loglik <- 0
  if (any(status == 1)) loglik <- loglik + sum(log(case_1_likelihood(V_0[as.integer(status) == 1],
                                                                     V_healthy[as.integer(status) == 1],
                                                                     T_obs[as.integer(status) == 1],
                                                                     a01, a02, a12, 
                                                                     A01, A02, A12)))
  
  if (any(as.integer(status) == 2)) loglik <- loglik + sum(log(case_2_likelihood(V_0[as.integer(status) == 2],
                                                                     V_healthy[as.integer(status) == 2],
                                                                     T_obs[as.integer(status) == 2],
                                                                     a01, a02, a12, 
                                                                     A01, A02, A12)))
  
  if (any(as.integer(status) == 3)) loglik <- loglik + sum(log(case_3_likelihood(V_0[as.integer(status) == 3],
                                                                     V_healthy[as.integer(status) == 3],
                                                                     V_ill[as.integer(status) == 3],
                                                                     T_obs[as.integer(status) == 3],
                                                                     a01, a02, a12, 
                                                                     A01, A02, A12)))
  
  if (any(as.integer(status) == 4)) loglik <- loglik + sum(log(case_4_likelihood(V_0[as.integer(status) == 4],
                                                                     V_healthy[as.integer(status) == 4],
                                                                     V_ill[as.integer(status) == 4],
                                                                     T_obs[as.integer(status) == 4],
                                                                     a01, a02, a12, 
                                                                     A01, A02, A12)))

  return(list(loglik = loglik,
              knots = knots, 
              degree = degree,
              hazards = list(a01 = a01,
                             a02 = a02,
                             a12 = a12,
                             A01 = A01,
                             A02 = A02,
                             A12 = A12)
              ))
}




do_likelihood_optim <- function(sim_dat, n_knots, degree, penalizer) {
  
  n_par = n_knots + degree - 1
  
  knots_a01 <- seq(min(sim_dat$V_0), max(sim_dat$T_obs), length.out = n_knots)
  knots_a02 <- knots_a01
  knots_a12 <- knots_a01
  
  knots = list(a01 = knots_a01,
               a02 = knots_a02,
               a12 = knots_a12)
  
  obj_fun <- function(x) {
    
    val <- -full_log_likehood(V_0 = sim_dat$V_0,
                              V_healthy = sim_dat$V_healthy,
                              V_ill = sim_dat$V_ill,
                              T_obs = sim_dat$T_obs,
                              status = sim_dat$status,
                              theta = list(a01 = x[1:n_par], 
                                           a02 = x[(n_par+1):(2*n_par)], 
                                           a12 = x[(2*n_par+1):(3*n_par)]),
                              degree, 
                              knots
    )$loglik - penalizer*(
      drop(t(x[1:n_par]) %*% pen_mat_m_splines(knots_a01) %*% x[1:n_par]) + 
        drop(t(x[(n_par+1):(2*n_par)]) %*% pen_mat_m_splines(knots_a02) %*% x[(n_par+1):(2*n_par)]) + 
        drop(t(x[(2*n_par+1):(3*n_par)]) %*% pen_mat_m_splines(knots_a12) %*% x[(2*n_par+1):(3*n_par)]))   
    
    if (!is.finite(val)) return(1e10)
    return(val)
  }
  
  x0 <- rep(1, 3*n_par)
  
  
  obj_fun(x0)
  
  
  res <- optim(
    par = x0,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = rep(0.0001, 3*n_par),
    upper = rep(13, 3*n_par)
  )
  
  res$par
  
  res_full <- full_log_likehood(V_0 = sim_dat$V_0,
                                V_healthy = sim_dat$V_healthy,
                                V_ill = sim_dat$V_ill,
                                T_obs = sim_dat$T_obs,
                                status = sim_dat$status,
                                theta = list(a01 = res$par[1:n_par], 
                                             a02 = res$par[(n_par+1):(2*n_par)], 
                                             a12 = res$par[(2*n_par+1):(3*n_par)]),
                                degree, 
                                knots)
  
  res_full
}
