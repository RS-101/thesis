plot_estimators <- function(res) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, 2))
  
  # Sort helper
  sort_xy <- function(x, y) {
    o <- order(x)
    list(x = x[o], y = y[o])
  }
  
  # F12
  xy <- sort_xy(res$s_eval, res$F12)
  plot(xy$x, xy$y, type = "s", col = "blue", lwd = 2,
       xlab = "s", ylab = expression(hat(F)[12](s)), main = "Estimator F12(s)")
  grid()
  
  # F13
  xy <- sort_xy(res$s_eval, res$F13)
  plot(xy$x, xy$y, type = "s", col = "red", lwd = 2,
       xlab = "s", ylab = expression(hat(F)[13](s)), main = "Estimator F13(s)")
  grid()
  
  # F
  xy <- sort_xy(res$s_eval, res$F)
  plot(xy$x, xy$y, type = "s", col = "purple", lwd = 2,
       xlab = "s", ylab = expression(hat(F)(s)), main = "Estimator F(s)")
  grid()
  
  # Hazards Λ12, Λ13, Λ23
  xy12 <- sort_xy(res$s_eval, res$Lambda12)
  xy13 <- sort_xy(res$s_eval, res$Lambda13)
  xy23 <- sort_xy(res$t_eval_for_23, res$Lambda23)
  
  plot(xy12$x, xy12$y, type = "s", col = "blue", lwd = 2,
       xlab = "s", ylab = "Cumulative hazards", main = "Λ estimators")
  lines(xy13$x, xy13$y, type = "s", col = "red", lwd = 2)
  lines(xy23$x, xy23$y, type = "s", col = "darkgreen", lwd = 2)
  legend("topleft",
         legend = c(expression(hat(Lambda)[12]), expression(hat(Lambda)[13]), expression(hat(Lambda)[23])),
         col = c("blue", "red", "darkgreen"), lwd = 2, bty = "n")
  grid()
}
