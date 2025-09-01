library(tidyverse)

source("likelihood_hazard_splines.R")
source("simulate_data.R")

sim_dat_all <- simulate_data(3000, censoring_time = 3)
sim_dat <- sim_dat_all$obs


do_likelihood_optim <- function(sim_dat, n_par, degree) {
  
  n_knots = n_par - degree  + 2 - 1
  
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
                              )$loglik
    
    if (!is.finite(val)) return(1e10)
    return(val)
  }
  
  x0 <- rep(1, 3*n_par)
  
  res <- optim(
    par = x0,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = rep(0.0001, 24),
    upper = rep(10, 24)
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

res_full <- do_likelihood_optim(sim_dat, n_par = 8, degree = 3)



# Collect your functions
funs <- list(
  a01 = res_full$hazards$a01,
  a02 = res_full$hazards$a02,
  a12 = res_full$hazards$a12
)

# Grid of x values
xvals <- seq(0, 3, length.out = 200)

# Evaluate each function
df <- lapply(names(funs), function(nm) {
  data.frame(
    x = xvals,
    y = sapply(xvals, funs[[nm]]),
    fun = nm
  )
}) %>% bind_rows()

# Plot
ggplot(df, aes(x = x, y = y, colour = fun)) +
  geom_line(size = 1) +
  labs(title = "Spline hazard",
       x = "x", y = "value") +
  theme_minimal()



