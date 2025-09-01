source("simulate_data.R")

dt <- simulate_illness_death_interval_censored_data(100)

k <- 4 # degree

g <- function(x) x^2

m <- 10 # number of knot
j <- 1 # knot index 

t_j <- 1 # knots

theta <- # spline parameters

sum(log(exp(-sum(g(theta)*I_j(data$L)))-exp(-sum(g(theta)*I_j(data$R))))) - kappa * 1 # int (sum(g(theta))*M''(j))

