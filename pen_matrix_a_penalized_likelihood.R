library(splines)

# Your existing B-spline penalty (second-derivative, cubic) -------------------
pen_mat_B <- function(inner_knots, ord = 4, derivs = 2) {
  # same construction pattern as your snippet
  knots <- sort(c(rep(range(inner_knots), ord - 1), inner_knots))
  d <- diff(inner_knots)                             # interval lengths
  mids <- inner_knots[-length(inner_knots)] + d / 2  # midpoints
  
  G_a <- splineDesign(knots, inner_knots, ord = ord, derivs = derivs)
  G_b <- G_a[-1, , drop = FALSE]
  G_a <- G_a[-nrow(G_a), , drop = FALSE]
  G_m <- splineDesign(knots, mids,        ord = ord, derivs = derivs)
  
  ( crossprod(d * G_a, G_a) +
      4 * crossprod(d * G_m, G_m) +
      crossprod(d * G_b, G_b) ) / 6
}

# M-spline penalty from the B-spline penalty ---------------------------------
pen_mat_M <- function(inner_knots, ord = 4, derivs = 2) {
  # 1) B-spline penalty for the same knots/order/derivative
  Omega_B <- pen_mat_B(inner_knots, ord = ord, derivs = derivs)
  
  # 2) Diagonal scaling D where M_i = alpha_i * N_i
  knots <- sort(c(rep(range(inner_knots), ord - 1), inner_knots))
  n <- length(knots) - ord
  alpha <- ord / (knots[(1:n) + ord] - knots[1:n])  # k / (t_{i+k} - t_i)
  D <- diag(alpha, nrow = n)
  
  # 3) Return Omega_M = D * Omega_B * D
  D %*% Omega_B %*% D
}
