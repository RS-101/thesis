library(qs)

source("joly/functions_estimation.R")
source("joly/functions_plot.R")
source("frydman/setup_cpp_from_data.R")
df_to_estimate <- list.files("experiments/data/")


for(i_df in seq_along(df_to_estimate)) {
  df <- df_to_estimate[i_df]
  full_df <- qread(paste0("experiments/data/", df))
  
  
  data_generation_plot <- full_df$true_data_generation$plot

  obs_data <- full_df$obs
  
  est_pl <- fit_idm(obs_data,
                    n_knots = 7, 
                    degree = 3, 
                    kappa_values = 10^(-2:3), 
                    verbose = T)
  
  est_npmle <- fit
  
  
}
