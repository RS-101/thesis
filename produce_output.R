library(qs)


files <- list.files("experiments/results",full.names = T)

for (i_df in seq_along(files)) {
  res <- qread(files[i_df])
  
  plot_estimators <- plot_cumhazards(res)
  
  plot_data <- res$full_df$true_data_generation$plot
  
}




