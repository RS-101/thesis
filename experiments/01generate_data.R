source("datageneration/functions_simulate_data.R")
library(qs)


data_size <- c(100, 500, 1000, 2000)

for (n in data_size) {
  sim_dat <- simulate_idm_constant_hazards(n = n, 
                                           a12 = 1/0.002,
                                           a23 = 1/0.005,
                                           a13 = 1/0.010)
  qsave(sim_dat, paste0("experiments/data/const_hazard_n_",n,".qs"))
  sim_dat <- simulate_idm_weibull(n,
                                  shape12 = 1.1, scale12 = 1/0.0008,
                                  shape13 = 1.8, scale13 = 1/0.0002,
                                  shape23 = 1.3, scale23 = 1/0.0016)
  qsave(sim_dat, paste0("experiments/data/wiebull_hazard_n_",n,".qs"))

  sim_dat <- simulate_idm_poly(n)
  qsave(sim_dat, paste0("experiments/data/general_poly_hazard_n_",n,".qs"))
}

rm(list = ls())
                                                                     