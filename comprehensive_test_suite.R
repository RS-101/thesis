library(Rcpp)
library(data.table)
library(ggplot2)
sourceCpp("optimizing_function_em.cpp")

source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")
source("frydman/functions_em.R")
source("frydman/helper_functions.R")

# ============================================================================
# COMPREHENSIVE TEST SUITE FOR Z_I SUM BUG
# ============================================================================

# Helper function to run a single test
run_single_test <- function(n, a12, a13, a23, seed, test_name, 
                            max_iter = 5, check_convergence = FALSE) {
  set.seed(seed)
  
  cat("\n")
  cat("========================================\n")
  cat("Test:", test_name, "\n")
  cat("  n =", n, ", a12 =", a12, ", a13 =", a13, ", a23 =", a23, "\n")
  cat("========================================\n")
  
  # Simulate data
  sim_data <- simulate_idm_constant_hazards(n = n, a12 = a12, a13 = a13, a23 = a23)
  
  # Convert to list format
  data_list <- setup_data_to_list_format(sim_data$obs, add_r_format = TRUE)
  
  # Quick stats
  cat("Data summary:\n")
  cat("  Total observations (N_star):", data_list$cpp_data$N_star, "\n")
  cat("  Case 3+4 (M):", data_list$cpp_data$M, "\n")
  cat("  Case 2 exact (K_tilde):", data_list$cpp_data$K_tilde, "\n")
  cat("  Case 2 interval (U):", data_list$cpp_data$U, "\n")
  cat("  Case 1 exact (J):", data_list$cpp_data$J, "\n")
  cat("  Case 1 interval (C):", data_list$cpp_data$C, "\n")
  cat("  Unique Q intervals (I):", data_list$cpp_data$I, "\n")
  cat("  Unique E_star (K):", data_list$cpp_data$K, "\n")
  cat("  Total parameters (I_mark):", data_list$cpp_data$I_mark, "\n")
  
  # Check c_k calculation
  cat("\nc_k verification:\n")
  cat("  sum(c_k) =", sum(data_list$cpp_data$c_k), "\n")
  cat("  K_tilde (expected) =", data_list$cpp_data$K_tilde, "\n")
  cat("  Match:", sum(data_list$cpp_data$c_k) == data_list$cpp_data$K_tilde, "\n")
  
  # Initialize
  set.seed(seed + 1000)
  z_init <- runif(data_list$cpp_data$I_mark)
  z_init <- z_init/sum(z_init)
  lambda_init <- runif(data_list$cpp_data$N, min = 0.1, max = 0.5)
  
  # Test C++ implementation
  cat("\nC++ Implementation (", max_iter, "iterations):\n", sep = "")
  tryCatch({
    res_cpp <- em_fit(
      make_model_data(data_list$cpp_data),
      z_init = z_init,
      lambda_init = lambda_init,
      max_iter = max_iter, 
      tol = 1e-6, 
      verbose = FALSE
    )
    
    cat("  sum(z_i) =", sum(res_cpp$z_i), "\n")
    cat("  max(z_i) =", max(res_cpp$z_i), "\n")
    cat("  min(z_i) =", min(res_cpp$z_i), "\n")
    cat("  sum(lambda_n) =", sum(res_cpp$lambda_n), "\n")
    cat("  Z sum test:", ifelse(abs(sum(res_cpp$z_i) - 1.0) < 1e-10, "PASS", "FAIL"), "\n")
    
    cpp_success <- TRUE
    cpp_z_sum <- sum(res_cpp$z_i)
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    cpp_success <<- FALSE
    cpp_z_sum <<- NA
  })
  
  # Test R implementation
  cat("\nR Implementation (", max_iter, "iterations):\n", sep = "")
  tryCatch({
    res_r <- do.call(
      em_estimate_raw,
      c(list(
        verbose = FALSE,
        max_iter = max_iter,
        tol = 1e-6,
        z_init = z_init,
        lambda_init = lambda_init
      ), data_list$r_data)
    )
    
    cat("  sum(z) =", sum(res_r$z), "\n")
    cat("  max(z) =", max(res_r$z), "\n")
    cat("  min(z) =", min(res_r$z), "\n")
    cat("  sum(lambda) =", sum(res_r$lambda), "\n")
    cat("  Z sum test:", ifelse(abs(sum(res_r$z) - 1.0) < 1e-10, "PASS", "FAIL"), "\n")
    
    r_success <- TRUE
    r_z_sum <- sum(res_r$z)
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    r_success <<- FALSE
    r_z_sum <<- NA
  })
  
  # Return test results
  list(
    test_name = test_name,
    n = n,
    a12 = a12, a13 = a13, a23 = a23,
    N_star = data_list$cpp_data$N_star,
    M = data_list$cpp_data$M,
    U = data_list$cpp_data$U,
    C = data_list$cpp_data$C,
    J = data_list$cpp_data$J,
    K_tilde = data_list$cpp_data$K_tilde,
    K = data_list$cpp_data$K,
    I = data_list$cpp_data$I,
    I_mark = data_list$cpp_data$I_mark,
    c_k_sum = sum(data_list$cpp_data$c_k),
    c_k_correct = sum(data_list$cpp_data$c_k) == data_list$cpp_data$K_tilde,
    cpp_success = cpp_success,
    cpp_z_sum = cpp_z_sum,
    cpp_z_ok = abs(cpp_z_sum - 1.0) < 1e-10,
    r_success = r_success,
    r_z_sum = r_z_sum,
    r_z_ok = abs(r_z_sum - 1.0) < 1e-10
  )
}

# ============================================================================
# TEST BATTERY
# ============================================================================

cat("\n")
cat("############################################################\n")
cat("# COMPREHENSIVE TEST SUITE\n")
cat("# Testing z_i normalization across various scenarios\n")
cat("############################################################\n")

results <- list()

# Test 1: Small sample, balanced hazards
results[[1]] <- run_single_test(
  n = 50, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016,
  seed = 123, test_name = "Small balanced"
)

# Test 2: Medium sample, balanced hazards
results[[2]] <- run_single_test(
  n = 250, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016,
  seed = 456, test_name = "Medium balanced"
)

# Test 3: Large sample, balanced hazards
results[[3]] <- run_single_test(
  n = 500, a12 = 0.001, a13 = 0.0003, a23 = 0.002,
  seed = 789, test_name = "Large balanced"
)

# Test 4: High illness rate (most go through state 2)
results[[4]] <- run_single_test(
  n = 200, a12 = 0.003, a13 = 0.0001, a23 = 0.002,
  seed = 111, test_name = "High illness rate"
)

# Test 5: High direct death rate (most skip state 2)
results[[5]] <- run_single_test(
  n = 200, a12 = 0.0001, a13 = 0.003, a23 = 0.002,
  seed = 222, test_name = "High direct death"
)

# Test 6: Very high death from illness (fast progression)
results[[6]] <- run_single_test(
  n = 200, a12 = 0.002, a13 = 0.0002, a23 = 0.01,
  seed = 333, test_name = "Fast progression through illness"
)

# Test 7: Very low death from illness (slow progression)
results[[7]] <- run_single_test(
  n = 200, a12 = 0.002, a13 = 0.0005, a23 = 0.0001,
  seed = 444, test_name = "Slow progression through illness"
)

# Test 8: Tiny sample (edge case)
results[[8]] <- run_single_test(
  n = 20, a12 = 0.001, a13 = 0.0003, a23 = 0.002,
  seed = 555, test_name = "Tiny sample"
)

# Test 9: Extra large sample
results[[9]] <- run_single_test(
  n = 1000, a12 = 0.0008, a13 = 0.0002, a23 = 0.0016,
  seed = 666, test_name = "Extra large sample"
)

# Test 10: Equal hazards
results[[10]] <- run_single_test(
  n = 150, a12 = 0.001, a13 = 0.001, a23 = 0.001,
  seed = 777, test_name = "Equal hazards"
)

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n")
cat("############################################################\n")
cat("# SUMMARY OF ALL TESTS\n")
cat("############################################################\n\n")

summary_df <- data.frame(
  Test = sapply(results, function(x) x$test_name),
  N = sapply(results, function(x) x$n),
  N_star = sapply(results, function(x) x$N_star),
  K_tilde = sapply(results, function(x) x$K_tilde),
  U = sapply(results, function(x) x$U),
  J = sapply(results, function(x) x$J),
  C = sapply(results, function(x) x$C),
  M = sapply(results, function(x) x$M),
  I_mark = sapply(results, function(x) x$I_mark),
  c_k_OK = sapply(results, function(x) ifelse(x$c_k_correct, "✓", "✗")),
  CPP_z_sum = sprintf("%.6f", sapply(results, function(x) x$cpp_z_sum)),
  CPP_OK = sapply(results, function(x) ifelse(x$cpp_z_ok, "✓", "✗")),
  R_z_sum = sprintf("%.6f", sapply(results, function(x) x$r_z_sum)),
  R_OK = sapply(results, function(x) ifelse(x$r_z_ok, "✓", "✗")),
  stringsAsFactors = FALSE
)

print(summary_df, row.names = FALSE)

# Check if all tests passed
all_cpp_ok <- all(sapply(results, function(x) x$cpp_z_ok))
all_r_ok <- all(sapply(results, function(x) x$r_z_ok))
all_c_k_ok <- all(sapply(results, function(x) x$c_k_correct))

cat("\n")
cat("Overall Results:\n")
cat("  c_k calculation: ", ifelse(all_c_k_ok, "ALL PASS ✓", "SOME FAIL ✗"), "\n")
cat("  C++ z_i sums:    ", ifelse(all_cpp_ok, "ALL PASS ✓", "SOME FAIL ✗"), "\n")
cat("  R z sums:        ", ifelse(all_r_ok, "ALL PASS ✓", "SOME FAIL ✗"), "\n")

if (all_cpp_ok && all_r_ok && all_c_k_ok) {
  cat("\n")
  cat("############################################################\n")
  cat("# SUCCESS! All tests passed! 🎉\n")
  cat("# The z_i normalization bug has been fixed.\n")
  cat("############################################################\n")
} else {
  cat("\n")
  cat("WARNING: Some tests failed. Further investigation needed.\n")
}
