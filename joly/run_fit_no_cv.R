# USING INTERCEPT = F IN SPLINES. IF WE WANT INTERCEPT WE NEED TO CHANGE SPLINE DESIGN FOR PEN MAT

source("datageneration/functions_simulate_data.R")
debugSource("joly/functions_likelihood.R")
debugSource("joly/functions_plot.R")
library(tidyverse)

simulated_data <- simulate_idm_weibull(1000, verbose = F)


res_full <- do_likelihood_optim(sim_dat = simulated_data$obs,
                                n_knots = 3, 
                                degree = 3,
                                penalizer = c(100, 100, 100))

plot(res_full$hazards, simulated_data$true_data_generation$hazards, cumulative = F, xlim = c(0,1))

