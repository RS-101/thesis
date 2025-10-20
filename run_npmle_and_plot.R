library(qs)
library(Rcpp)

source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")
source("experiments/plot_function.R")

sim_data <- simulate_idm_constant_hazards(n = 250, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016)


est_npmle <- get_npmle(sim_data$obs, 
                       max_iter = 1000, tol = 0.0001, verbose = TRUE)

plot_estimators <- plot_cumhazards(list(full_df = sim_data, npmle = est_npmle))

est_npmle$raw_em_res$z_i %>% sum()

plot_estimators_gg(est_npmle$estimators_point)
plot_estimators$plot
