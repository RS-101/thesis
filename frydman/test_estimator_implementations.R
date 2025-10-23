# Test script to compare two implementations of estimators
# Implementation 1: calc_F_and_hazards() from helper_functions.R
# Implementation 2: Individual calc_* functions from estimation_of_A_to_test_against.R

library(testthat)
library(qs)
library(Rcpp)

# Source the implementations
source("frydman/helper_functions.R")
source("frydman/estimation_of_A_functions_only.R")
source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")

cat("=== Testing Estimator Implementations ===\n\n")

# Generate test data
set.seed(123)
cat("Generating test data...\n")
sim_data <- simulate_idm_constant_hazards(n = 100, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016)
data_list <- setup_data_to_list_format(sim_data$obs)

# Run NPMLE estimation
cat("Running NPMLE estimation...\n")
est_npmle <- get_npmle(NULL, data_list, max_iter = 1000, tol = 0.01, verbose = FALSE)

# Extract results
z_i_vec <- est_npmle$raw_em_res$z_i
lambda_n <- est_npmle$raw_em_res$lambda
Q_i <- data_list$Q_i
Q_i_mark <- data_list$Q_i_mark
t_star_n <- data_list$t_star_n
I <- data_list$I
I_mark <- data_list$I_mark

cat("\nData dimensions:\n")
cat("  I =", I, "\n")
cat("  I_mark =", I_mark, "\n")
cat("  N (lambda_n length) =", length(lambda_n), "\n")
cat("  z_i length =", length(z_i_vec), "\n")

# Define test grid
max_time <- max(Q_i[is.finite(Q_i[,2]), 2], Q_i_mark[is.finite(Q_i_mark)], na.rm = TRUE)
grid_points <- seq(0, max_time, length.out = 100)

cat("\n=== Implementation 1: calc_F_and_hazards() ===\n")
impl1 <- calc_F_and_hazards(
  grid_points = grid_points,
  z_i = z_i_vec,
  lambda_n = lambda_n,
  Q_i = Q_i,
  Q_i_mark = Q_i_mark,
  t_star_n = t_star_n,
  E_star = data_list$E_star
)

cat("Results computed successfully\n")
cat("  F12 range: [", min(impl1$as_points$F12), ",", max(impl1$as_points$F12), "]\n")
cat("  F13 range: [", min(impl1$as_points$F13), ",", max(impl1$as_points$F13), "]\n")
cat("  F range: [", min(impl1$as_points$F), ",", max(impl1$as_points$F), "]\n")

cat("\n=== Implementation 2: Individual calc_* functions ===\n")
F12_impl2 <- calc_F_12(z_i_vec, Q_i, I)$fill
F13_impl2 <- calc_F_13(z_i_vec, Q_i_mark, I, I_mark)
F_impl2 <- calc_F(z_i_vec, Q_i, Q_i_mark, I, I_mark)
A12_impl2 <- calc_A_12(z_i_vec, Q_i, Q_i_mark, I, I_mark)
A13_impl2 <- calc_A_13(z_i_vec, Q_i, Q_i_mark, I, I_mark)
A23_impl2 <- calc_A_23(lambda_n, t_star_n)

# Evaluate at grid points
F12_impl2_vals <- sapply(grid_points, F12_impl2)
F13_impl2_vals <- sapply(grid_points, F13_impl2)
F_impl2_vals <- sapply(grid_points, F_impl2)
A12_impl2_vals <- sapply(grid_points, A12_impl2)
A13_impl2_vals <- sapply(grid_points, A13_impl2)
A23_impl2_vals <- sapply(grid_points, A23_impl2)

cat("Results computed successfully\n")
cat("  F12 range: [", min(F12_impl2_vals), ",", max(F12_impl2_vals), "]\n")
cat("  F13 range: [", min(F13_impl2_vals), ",", max(F13_impl2_vals), "]\n")
cat("  F range: [", min(F_impl2_vals), ",", max(F_impl2_vals), "]\n")

cat("\n=== Comparison of Results ===\n\n")

# Compare F12
diff_F12 <- impl1$as_points$F12 - F12_impl2_vals
cat("F12:\n")
cat("  Max absolute difference:", max(abs(diff_F12)), "\n")
cat("  Mean absolute difference:", mean(abs(diff_F12)), "\n")
cat("  Relative error (max):", max(abs(diff_F12) / (abs(F12_impl2_vals) + 1e-10)), "\n")

# Compare F13
diff_F13 <- impl1$as_points$F13 - F13_impl2_vals
cat("\nF13:\n")
cat("  Max absolute difference:", max(abs(diff_F13)), "\n")
cat("  Mean absolute difference:", mean(abs(diff_F13)), "\n")
cat("  Relative error (max):", max(abs(diff_F13) / (abs(F13_impl2_vals) + 1e-10)), "\n")

# Compare F (total)
diff_F <- impl1$as_points$F - F_impl2_vals
cat("\nF (total):\n")
cat("  Max absolute difference:", max(abs(diff_F)), "\n")
cat("  Mean absolute difference:", mean(abs(diff_F)), "\n")
cat("  Relative error (max):", max(abs(diff_F) / (abs(F_impl2_vals) + 1e-10)), "\n")

# Compare Lambda12
diff_A12 <- impl1$as_points$Lambda12 - A12_impl2_vals
cat("\nLambda12:\n")
cat("  Max absolute difference:", max(abs(diff_A12)), "\n")
cat("  Mean absolute difference:", mean(abs(diff_A12)), "\n")
cat("  Relative error (max):", max(abs(diff_A12) / (abs(A12_impl2_vals) + 1e-10)), "\n")

# Compare Lambda13
diff_A13 <- impl1$as_points$Lambda13 - A13_impl2_vals
cat("\nLambda13:\n")
cat("  Max absolute difference:", max(abs(diff_A13)), "\n")
cat("  Mean absolute difference:", mean(abs(diff_A13)), "\n")
cat("  Relative error (max):", max(abs(diff_A13) / (abs(A13_impl2_vals) + 1e-10)), "\n")

# Compare Lambda23
diff_A23 <- impl1$as_points$Lambda23 - A23_impl2_vals
cat("\nLambda23:\n")
cat("  Max absolute difference:", max(abs(diff_A23)), "\n")
cat("  Mean absolute difference:", mean(abs(diff_A23)), "\n")
cat("  Relative error (max):", max(abs(diff_A23) / (abs(A23_impl2_vals) + 1e-10)), "\n")

cat("\n=== Summary ===\n")
tolerance <- 1e-6
all_close <- all(
  abs(diff_F12) < tolerance,
  abs(diff_F13) < tolerance,
  abs(diff_F) < tolerance,
  abs(diff_A12) < tolerance,
  abs(diff_A13) < tolerance,
  abs(diff_A23) < tolerance
)

if (all_close) {
  cat("✓ All implementations match within tolerance (", tolerance, ")\n", sep = "")
} else {
  cat("✗ Implementations differ beyond tolerance (", tolerance, ")\n", sep = "")
  cat("\nLargest discrepancies:\n")
  max_diffs <- c(
    F12 = max(abs(diff_F12)),
    F13 = max(abs(diff_F13)),
    F = max(abs(diff_F)),
    Lambda12 = max(abs(diff_A12)),
    Lambda13 = max(abs(diff_A13)),
    Lambda23 = max(abs(diff_A23))
  )
  print(sort(max_diffs, decreasing = TRUE))
}

cat("\n=== Visual Comparison (First 10 points) ===\n")
comparison_df <- data.frame(
  grid = grid_points[1:10],
  F12_impl1 = impl1$as_points$F12[1:10],
  F12_impl2 = F12_impl2_vals[1:10],
  F13_impl1 = impl1$as_points$F13[1:10],
  F13_impl2 = F13_impl2_vals[1:10],
  A12_impl1 = impl1$as_points$Lambda12[1:10],
  A12_impl2 = A12_impl2_vals[1:10]
)
print(comparison_df)

cat("\n=== Test Complete ===\n")
