# Final verification: Demonstrate the bug and proposed fix

library(Rcpp)
library(qs)

source("frydman/helper_functions.R")
source("frydman/estimation_of_A_functions_only.R")
source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")

cat("=== FINAL VERIFICATION OF BUG AND FIX ===\n\n")

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

max_time <- max(Q_i[is.finite(Q_i[,2]), 2], na.rm = TRUE)
grid_points <- seq(0, max_time, length.out = 100)

cat("1. Current Implementation (BUGGY)\n")
cat("   --------------------------------\n")
impl_buggy <- calc_F_and_hazards(
  grid_points = grid_points,
  z_i = z_i_vec,
  lambda_n = lambda_n,
  Q_i = Q_i,
  Q_i_mark = Q_i_mark,
  t_star_n = t_star_n,
  E_star = data_list$E_star
)

cat("2. Reference Implementation (CORRECT)\n")
cat("   -----------------------------------\n")
A12_correct <- calc_A_12(z_i_vec, Q_i, Q_i_mark, I, I_mark)
A13_correct <- calc_A_13(z_i_vec, Q_i, Q_i_mark, I, I_mark)
A12_correct_vals <- sapply(grid_points, A12_correct)
A13_correct_vals <- sapply(grid_points, A13_correct)

cat("3. Proposed Fix (CORRECTED Implementation 1)\n")
cat("   ------------------------------------------\n")

# Manually implement the corrected version
step_cdf <- function(grid_points, times, masses, ...) {
  time_order <- order(times)
  times <- times[time_order]
  masses <- masses[time_order]
  cs <- pmax(cumsum(masses),0)
  idx <- findInterval(grid_points, times, ...)
  c(0,cs)[idx+1]
}

# Corrected Lambda12 calculation
F12_at_l <- step_cdf(Q_i[,1]-1e-6, times = Q_i[1:I, 2], masses = z_i_vec[1:I])
F13_at_l <- step_cdf(Q_i[,1]-1e-6, times = Q_i_mark, masses = z_i_vec[(I+1):I_mark])
F_total_at_l <- F12_at_l + F13_at_l  # KEY FIX: Use F_total, not just F12

denom12_fixed <- 1 - F_total_at_l
term12_fixed <- ifelse(denom12_fixed > 0, z_i_vec[1:I] / denom12_fixed, 0)
A12_fixed <- step_cdf(grid_points, times = Q_i[,2], term12_fixed)

# Corrected Lambda13 calculation
F12_at_e <- step_cdf(Q_i_mark-1e-6, times = Q_i[1:I, 2], masses = z_i_vec[1:I])
F13_at_e <- step_cdf(Q_i_mark-1e-6, times = Q_i_mark, masses = z_i_vec[(I+1):I_mark])
F_total_at_e <- F12_at_e + F13_at_e  # KEY FIX: Use F_total, not just F13

denom13_fixed <- 1 - F_total_at_e
term13_fixed <- ifelse(denom13_fixed > 0, z_i_vec[(I+1):I_mark] / denom13_fixed, 0)
A13_fixed <- step_cdf(grid_points, times = Q_i_mark, term13_fixed)

cat("\n=== COMPARISON RESULTS ===\n\n")

# Compare Lambda12
diff_buggy_vs_correct <- impl_buggy$as_points$Lambda12 - A12_correct_vals
diff_fixed_vs_correct <- A12_fixed - A12_correct_vals

cat("Lambda12:\n")
cat("  Buggy vs Correct - Max diff:", max(abs(diff_buggy_vs_correct)), "\n")
cat("  Fixed vs Correct - Max diff:", max(abs(diff_fixed_vs_correct)), "\n")
cat("  Improvement:", max(abs(diff_buggy_vs_correct)) - max(abs(diff_fixed_vs_correct)), "\n\n")

# Compare Lambda13
diff_buggy_vs_correct_13 <- impl_buggy$as_points$Lambda13 - A13_correct_vals
diff_fixed_vs_correct_13 <- A13_fixed - A13_correct_vals

cat("Lambda13:\n")
cat("  Buggy vs Correct - Max diff:", max(abs(diff_buggy_vs_correct_13)), "\n")
cat("  Fixed vs Correct - Max diff:", max(abs(diff_fixed_vs_correct_13)), "\n")
cat("  Improvement:", max(abs(diff_buggy_vs_correct_13)) - max(abs(diff_fixed_vs_correct_13)), "\n\n")

cat("=== VERIFICATION ===\n\n")

tolerance <- 1e-10
if (max(abs(diff_fixed_vs_correct)) < tolerance && 
    max(abs(diff_fixed_vs_correct_13)) < tolerance) {
  cat("✅ SUCCESS! The fix makes Implementation 1 match Implementation 2 perfectly.\n\n")
  cat("The bug was confirmed and the fix is validated.\n\n")
} else {
  cat("⚠️  The fix reduces the error but doesn't eliminate it completely.\n")
  cat("   Further investigation may be needed.\n\n")
}

cat("=== SUMMARY ===\n\n")
cat("BUG IDENTIFIED:\n")
cat("  • calc_F_and_hazards() uses F12 in Lambda12 denominator (should use F_total)\n")
cat("  • calc_F_and_hazards() uses F13 in Lambda13 denominator (should use F_total)\n\n")
cat("FIX REQUIRED:\n")
cat("  • Change: denom12 <- 1 - F12_at_l_i\n")
cat("  • To:     denom12 <- 1 - (F12_at_l_i + F13_at_l_i)\n\n")
cat("  • Change: denom13 <- 1 - F_total_at_e_k\n")
cat("  • To:     denom13 <- 1 - (F12_at_e_k + F13_at_e_k)\n\n")
cat("THEORETICAL JUSTIFICATION:\n")
cat("  The risk set = {individuals who haven't experienced ANY event yet}\n")
cat("  = 1 - P(any event by time t) = 1 - [F12(t) + F13(t)] = 1 - F_total(t)\n\n")

cat("=== Test Complete ===\n")
