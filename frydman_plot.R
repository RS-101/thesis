# ggplot2 rewrite (uses ggplot2 + patchwork)
library(ggplot2)
library(dplyr)
library(patchwork)

plot_estimators_gg <- function(res) {
  # Common theme
  thm <- theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_line(),
          plot.title = element_text(face = "bold"))
  
  # F12
  df_F12 <- data.frame(s = res$grid_points, y = res$F12) |> arrange(s)
  p1 <- ggplot(df_F12, aes(x = s, y = y)) +
    geom_step(linewidth = 1, color = "blue") +
    labs(x = "s", y = expression(hat(F)[12](s)), title = "Estimator F12(s)") +
    thm
  
  # F13
  df_F13 <- data.frame(s = res$grid_points, y = res$F13) |> arrange(s)
  p2 <- ggplot(df_F13, aes(x = s, y = y)) +
    geom_step(linewidth = 1, color = "red") +
    labs(x = "s", y = expression(hat(F)[13](s)), title = "Estimator F13(s)") +
    thm
  
  # F
  df_F <- data.frame(s = res$grid_points, y = res$F) |> arrange(s)
  p3 <- ggplot(df_F, aes(x = s, y = y)) +
    geom_step(linewidth = 1, color = "purple") +
    labs(x = "s", y = expression(hat(F)(s)), title = "Estimator F(s)") +
    thm
  
  # Hazards Λ12, Λ13, Λ23
  df_haz <- dplyr::bind_rows(
    data.frame(s = res$grid_points,         y = res$Lambda12, which = "L12"),
    data.frame(s = res$grid_points,         y = res$Lambda13, which = "L13"),
    data.frame(s = res$grid_points,         y = res$Lambda23, which = "L23")
  ) |> arrange(s)
  
  p4 <- ggplot(df_haz, aes(x = s, y = y, color = which)) +
    geom_step(linewidth = 1, alpha = 0.5) +
    scale_color_manual(
      values = c(L12 = "blue", L13 = "red", L23 = "darkgreen"),
      labels = c(
        L12 = expression(hat(Lambda)[12]),
        L13 = expression(hat(Lambda)[13]),
        L23 = expression(hat(Lambda)[23])
      )
    ) +
    labs(x = "s", y = "Cumulative hazards", color = NULL, title = "Λ estimators") +
    thm
  
  # Arrange 2x2 grid
  (p1 | p2) / (p3 | p4)
}
