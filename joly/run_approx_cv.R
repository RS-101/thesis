source("joly/functions_likelihood.R")
source("datageneration/functions_simulate_data.R")
source("joly/functions_plot.R")

dat <- simulate_idm_constant_hazards(100)
res_cv <- fit_idm(dat$obs)


res_greedy <- fit_idm_greedy(dat$obs)



res_cv$fit$kappa_term


res_cv$fit$knots

plot(res_cv$fit$hazards, true_hazards = dat$true_data_generation$hazards)


# res <- do_likelihood_optim(dat$obs, kappa_term = c(1,10,0.1), n_knots = 7)


