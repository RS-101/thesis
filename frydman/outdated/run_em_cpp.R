library(Rcpp)
source("datageneration/functions_simulate_data.R")
sourceCpp("frydman/functions_em.cpp")
source("frydman/setup_cpp_from_data.R")
source("frydman/functions_extract_estimators.R")
source("frydman/functions_plot.R")
source("frydman/setup_for_R_from_data.R")

set.seed(100)
sim_dat <- simulate_idm_weibull(2000,
                                shape12 = 1.1, scale12 = 1/0.0008,
                                shape13 = 1.8, scale13 = 1/0.0002,
                                shape23 = 1.3, scale23 = 1/0.0016)

sim_dat$true_data_generation$plot
