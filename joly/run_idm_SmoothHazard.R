library(SmoothHazard)
library(prodlim)
d_full <- simulate_idm_weibull(1000)

d <- d_full$obs
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
  CV        = T,          # tune smoothing (kappa) via approx. CV
  conf.int  = TRUE
)

print(fit_spl)

plot(fit_spl)

h_weibull <- function(t, shape, scale) {
  t <- pmax(as.numeric(t), .Machine$double.eps)
  (shape / scale) * (t / scale)^(shape - 1)
}

a12 <- function(t) h_weibull(t, fit_spl$modelPar[1], fit_spl$modelPar[2])  # 1 -> 2
a13 <- function(t) h_weibull(t, fit_spl$modelPar[3], fit_spl$modelPar[4])  # 1 -> 3
a23 <- function(t) h_weibull(t, fit_spl$modelPar[5], fit_spl$modelPar[6])  # 2 -> 3 

funs <- list(
  a01 = a12,
  a02 = a13,
  a12 = a23,
  true_a01 = d_full$true_data_generation$hazards$a12,
  true_a02 = d_full$true_data_generation$hazards$a13,
  true_a12 = d_full$true_data_generation$hazards$a23
)

# Grid of x values
xvals <- seq(0, 1, length.out = 200)



# Evaluate each function
df <- lapply(names(funs), function(nm) {
  data.frame(
    x = xvals,
    y = sapply(xvals, funs[[nm]]),
    fun = nm
  )
}) %>% bind_rows()

df <- df %>% mutate(type = ifelse(str_detect(fun, "true"), "true", "estimate"))


# Plot
ggplot(df, aes(x = x, y = y, colour = fun)) +
  geom_line(size = 1, aes(linetype = type)) +
  labs(title = "Spline hazard",
       x = "x", y = "value") +
  theme_minimal()


