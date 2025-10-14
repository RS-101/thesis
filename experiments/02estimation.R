library(qs)
library(Rcpp)


source("joly/functions_estimation.R")
source("joly/functions_plot.R")
source("frydman/wrapper.R")


df_to_estimate <- list.files("experiments/data/")[1]


for(i_df in seq_along(df_to_estimate)) {
  df <- df_to_estimate[i_df]
  full_df <- qread(paste0("experiments/data/", df))
  
  
  data_generation_plot <- full_df$true_data_generation$plot

  obs_data <- full_df$obs
  
  est_npmle <- get_npmle(obs_data, max_iter = 100, tol = 0.001, verbose = T)
  
  est_pl <- fit_idm(obs_data,
                    n_knots = 7, 
                    degree = 3, 
                    kappa_values = c(0.1,10,100), 
                    verbose = T)
  
  qsave(list(npmle = est_npmle,
             penmle = est_pl, 
             full_df),
        file = paste0("experiments/results/estimate_",df))
}
