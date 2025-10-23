# Deep dive into Lambda12 differences
library(ggplot2)
library(Rcpp)
library(qs)

source("frydman/helper_functions.R")
source("frydman/estimation_of_A_functions_only.R")
source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")

cat("=== Investigating Lambda12 Differences ===\n\n")

# Generate test data
set.seed(123)
sim_data <- simulate_idm_constant_hazards(n = 100, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016)
data_list <- setup_data_to_list_format(sim_data$obs)
est_npmle <- get_npmle(NULL, data_list, max_iter = 1000, tol = 0.01, verbose = FALSE)

z_i_vec <- est_npmle$raw_em_res$z_i
lambda_n <- est_npmle$raw_em_res$lambda
Q_i <- data_list$Q_i
Q_i_mark <- data_list$Q_i_mark
t_star_n <- data_list$t_star_n
I <- data_list$I
I_mark <- data_list$I_mark

cat("Key parameters:\n")
cat("  I =", I, "\n")
cat("  Sum of z_i[1:I] =", sum(z_i_vec[1:I]), "\n")
cat("  Number of unique jump points in Q_i =", length(unique(Q_i[,2])), "\n\n")

# Look at the F12 calculation in detail
max_time <- max(Q_i[is.finite(Q_i[,2]), 2], na.rm = TRUE)
test_points <- unique(sort(c(
  seq(0, max_time, length.out = 50),
  Q_i[,2][Q_i[,2] < Inf]
)))

cat("=== Testing F12 at jump points ===\n\n")

# Implementation 1
impl1 <- calc_F_and_hazards(
  grid_points = test_points,
  z_i = z_i_vec,
  lambda_n = lambda_n,
  Q_i = Q_i,
  Q_i_mark = Q_i_mark,
  t_star_n = t_star_n,
  E_star = data_list$E_star
)

# Implementation 2
F_impl2 <- calc_F(z_i_vec, Q_i, Q_i_mark, I, I_mark)
A12_impl2 <- calc_A_12(z_i_vec, Q_i, Q_i_mark, I, I_mark)

# Evaluate
F_impl2_vals <- sapply(test_points, F_impl2)
A12_impl2_vals <- sapply(test_points, A12_impl2)

# Check F values at l_i (left endpoints)
cat("Checking F at left endpoints (l_i - epsilon):\n")
F12_at_l_impl1 <- impl1$as_functions$F12(Q_i[,1] - 1e-6)
F_at_l_impl2 <- F_impl2(Q_i[,1] - 1e-10)
F12_at_l_impl2 <- sapply(Q_i[,1] - 1e-10, function(s) {
  sum(z_i_vec[1:I] * as.numeric(Q_i[1:I,2] <= s))
})

cat("  First 5 comparisons:\n")
comparison <- data.frame(
  i = 1:5,
  l_i = Q_i[1:5, 1],
  F_impl1 = F12_at_l_impl1[1:5],
  F_impl2 = F12_at_l_impl2[1:5],
  diff = F12_at_l_impl1[1:5] - F12_at_l_impl2[1:5]
)
print(comparison)

cat("\n=== Examining Lambda12 calculation ===\n\n")

# Let's manually compute Lambda12 using both methods for a few points
cat("Lambda12 formula: sum_{i: r_i <= s} z_i / (1 - F(l_i-))\n\n")

# Pick a test point after several jumps
s_test <- Q_i[10, 2]
cat("Testing at s =", s_test, "(r_10)\n\n")

# Implementation 1 approach
F12_at_l_i_impl1 <- sapply(Q_i[,1], impl1$as_functions$F12)
z_i_12 <- z_i_vec[1:I]
idx <- Q_i[,2] <= s_test
denom_impl1 <- 1 - F12_at_l_i_impl1
term_impl1 <- ifelse(denom_impl1 > 0, z_i_12 / denom_impl1, 0)
Lambda12_impl1_manual <- sum(term_impl1[idx])

cat("Implementation 1 (manual calc):\n")
cat("  Indices where r_i <= s:", which(idx), "\n")
cat("  Lambda12 =", Lambda12_impl1_manual, "\n\n")

# Implementation 2 approach  
F_hat_at_l <- F_impl2(Q_i[,1] - 1e-10)
denom_impl2 <- 1 - F_hat_at_l
term_impl2 <- ifelse(denom_impl2 > 0, z_i_12 / denom_impl2, 0)
Lambda12_impl2_manual <- sum(term_impl2[idx])

cat("Implementation 2 (manual calc):\n")
cat("  F(l_i-) for first few i:", F_hat_at_l[1:5], "\n")
cat("  Lambda12 =", Lambda12_impl2_manual, "\n\n")

cat("Difference:", Lambda12_impl1_manual - Lambda12_impl2_manual, "\n\n")

# Now let's check if the issue is in F12 vs F (total)
cat("=== Key insight: Implementation 2 uses F_total, not F12 ===\n\n")

# Implementation 2 uses F (total) = F12 + F13 in the denominator
# Let's check if this is the source of difference

F_total_at_l <- sapply(Q_i[,1] - 1e-6, impl1$as_functions$F)
F12_at_l <- sapply(Q_i[,1] - 1e-6, impl1$as_functions$F12)

cat("Comparison of denominators (first 5):\n")
denom_comparison <- data.frame(
  i = 1:5,
  using_F_total = 1 - F_total_at_l[1:5],
  using_F12 = 1 - F12_at_l[1:5],
  difference = F_total_at_l[1:5] - F12_at_l[1:5]
)
print(denom_comparison)

cat("\n=== DIAGNOSIS ===\n\n")
cat("The difference is because:\n")
cat("  • Implementation 1 (calc_F_and_hazards) uses F12 only: F12(l_i-)\n")
cat("  • Implementation 2 (calc_A_12) uses F total: F(l_i-) = F12(l_i-) + F13(l_i-)\n\n")
cat("This is a CONCEPTUAL DIFFERENCE, not a numerical precision issue!\n\n")

cat("According to the theory, Lambda12 should use:\n")
cat("  The risk set = 1 - F_total(l_i-), since anyone who has experienced\n")
cat("  either event 2 or event 3 before l_i is not at risk.\n\n")

cat("Therefore: Implementation 2 (using F total) is THEORETICALLY CORRECT\n")
cat("           Implementation 1 (using F12 only) is INCORRECT\n\n")

# Verify by checking against theoretical formulation
cat("=== Recommendation ===\n")
cat("Fix calc_F_and_hazards() to use F_total instead of F12 in Lambda12 calculation.\n")
cat("Similarly, Lambda13 should also use F_total.\n")
