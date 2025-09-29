# USING INTERCEPT = F IN SPLINES. IF WE WANT INTERCEPT WE NEED TO CHANGE SPLINE DESIGN FOR PEN MAT


library(tidyverse)

debugSource("likelihood_hazard_splines.R")



res_full <- do_likelihood_optim(res_w_ic$obs, n_knots = 6, degree = 3, penalizer = 0.5)



# Collect your functions
funs <- list(
  a01 = res_full$hazards$a01,
  a02 = res_full$hazards$a02,
  a12 = res_full$hazards$a12,
  true_a01 = res$hazards$a12,
  true_a02 = res$hazards$a13,
  true_a12 = res$hazards$a23
)

# Grid of x values
xvals <- seq(0, 13, length.out = 200)



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
  theme_minimal() + 
  ylim(0,10)

