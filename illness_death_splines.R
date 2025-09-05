# USING INTERCEPT = F IN SPLINES. IF WE WANT INTERCEPT WE NEED TO CHANGE SPLINE DESIGN FOR PEN MAT


library(tidyverse)



debugSource("likelihood_hazard_splines.R")
debugSource("pen_matrix_a_penalized_likelihood.R")
source("simulate_data.R")

sim_dat_all <- simulate_data(1000, censoring_time = 5)
sim_dat <- sim_dat_all$data$obs

sim_dat %>% summary()

do_likelihood_optim <- function(sim_dat, n_par, degree) {
  
  n_knots = n_par - degree  + 2
  
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
                              )$loglik - 0.1*(
      drop(t(x[1:n_par]) %*% pen_mat_M(knots_a01[2:n_knots]) %*% x[1:n_par]) - 
      drop(t(x[(n_par+1):(2*n_par)]) %*% pen_mat_M(knots_a02[2:n_knots]) %*% x[(n_par+1):(2*n_par)]) - 
      drop(t(x[(2*n_par+1):(3*n_par)]) %*% pen_mat_M(knots_a12[2:n_knots]) %*% x[(2*n_par+1):(3*n_par)]))   
                              
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
    upper = rep(10, 3*n_par)
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

res_full <- do_likelihood_optim(sim_dat, n_par = 7, degree = 3)



# Collect your functions
funs <- list(
  a01 = res_full$hazards$a01,
  a02 = res_full$hazards$a02,
  a12 = res_full$hazards$a12,
  true_a01 = sim_dat_all$hazards$a01,
  true_a02 = sim_dat_all$hazards$a02,
  true_a12 = sim_dat_all$hazards$a12
)

# Grid of x values
xvals <- seq(0, 5, length.out = 200)



# Evaluate each function
df <- lapply(names(funs), function(nm) {
  data.frame(
    x = xvals,
    y = sapply(xvals, funs[[nm]]),
    fun = nm
  )
}) %>% bind_rows()

df <- df %>% mutate(type = ifelse(str_detect(fun, "true"), "true", "estimate"))


# Plot
ggplot(df, aes(x = x, y = y, colour = fun)) +
  geom_line(size = 1, aes(linetype = type)) +
  labs(title = "Spline hazard",
       x = "x", y = "value") +
  theme_minimal() 



