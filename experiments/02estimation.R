library(qs)

df_to_estimate <- list.files("experiments/data/")


for(i_df in seq_along(df_to_estimate)) {
  df <- df_to_estimate[i_df]
  full_df <- qread(paste0("experiments/data/", df))
  
  
  data_generation_plot <- full_df$true_data_generation$plot

  obs_data <- full_df$obs
  
  
  
}
