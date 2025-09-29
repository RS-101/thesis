library(SmoothHazard)
library(prodlim)

d <- res_w_ic$obs

# Build interval-censored illness (0->1)
d$seen_ill  <- !is.na(d$V_ill)
d$L <- d$V_healthy
d$R <- ifelse(d$seen_ill, d$V_ill, d$T_obs)

# Death indicator used for 0->2 and 1->2
d$seen_exit <- grepl("died", d$status)

# Illness–death model with M-splines (penalized likelihood)
fit_spl <- idm(
  formula01 = Hist(time = list(L, R), event = seen_ill) ~ 1,
  formula02 = Hist(time = T_obs,       event = seen_exit) ~ 1,
  formula12 = Hist(time = T_obs,       event = seen_exit) ~ 1,
  data      = d,
  method    = "Splines",     # spline baseline hazards
  n.knots   = c(9, 9, 9),    # reasonable smoothness
  CV        = F,          # tune smoothing (kappa) via approx. CV
  conf.int  = TRUE
)

print(fit_spl)

# Plot baseline transition intensities
op <- par(mfrow = c(1, 3), mar = c(4,4,2,1))
plot(fit_spl, transition = "01", conf.int = TRUE, xlab = "Time", ylab = "h01(t)")
title("0 -> 1 (illness)")
plot(fit_spl, transition = "02", conf.int = TRUE, xlab = "Time", ylab = "h02(t)")
title("0 -> 2 (death without illness)")
plot(fit_spl, transition = "12", conf.int = TRUE, xlab = "Time", ylab = "h12(t)")
title("1 -> 2 (death after illness)")
par(op)
