
list_1 <- setup_frydman(data)

res_em <- do.call(
  em_estimate,
  c(list(
    verbose = TRUE,
    max_iter = 300
  ),list_1)
)


estimators <- calc_F_and_hazards(
  grid_points = seq(0, 3, by = 0.01),
  z_i = res_em$z, lambda = res_em$lambda, 
  list_1$Q_full, list_1$T_star, list_1$E_star
)

plot_estimators_gg(estimators)



