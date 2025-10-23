library(Rcpp)
library(data.table)
sourceCpp("optimizing_function_em.cpp")

source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")
source("frydman/functions_em.R")
source("frydman/helper_functions.R")

# ============================================================================
# DIAGNOSTIC SCRIPT: Find the bug causing z_i to not sum to 1
# ============================================================================

# Create a comprehensive testing function
test_z_sum <- function(sim_data, test_name = "Test") {
  cat("\n")
  cat("==========================================================\n")
  cat("Testing:", test_name, "\n")
  cat("==========================================================\n")
  
  data_list <- setup_data_to_list_format(sim_data$obs, add_r_format = TRUE)
  
  # Print diagnostic info about the data
  cat("\nData Structure:\n")
  cat("  N (unique t_star_n):", data_list$cpp_data$N, "\n")
  cat("  M (case 3+4):", data_list$cpp_data$M, "\n")
  cat("  U (case 2, non-exact):", data_list$cpp_data$U, "\n")
  cat("  C (case 1, non-exact):", data_list$cpp_data$C, "\n")
  cat("  J (case 1, exact):", data_list$cpp_data$J, "\n")
  cat("  K_tilde (case 2, exact):", data_list$cpp_data$K_tilde, "\n")
  cat("  K (unique E_star):", data_list$cpp_data$K, "\n")
  cat("  I (Q_i intervals):", data_list$cpp_data$I, "\n")
  cat("  I_mark (I + K):", data_list$cpp_data$I_mark, "\n")
  cat("  N_star (total count):", data_list$cpp_data$N_star, "\n")
  cat("  W (M + U):", data_list$cpp_data$W, "\n")
  cat("  M_mark (M + U + C):", data_list$cpp_data$M_mark, "\n")
  
  # Check c_k calculation
  cat("\n")
  cat("c_k Diagnostic:\n")
  e_k <- data_list$cpp_data$e_k
  t_u <- data_list$cpp_data$t_u
  E_star <- data_list$cpp_data$E_star
  
  # Method 1: count only e_k (used in functions_em.R)
  c_k_method1 <- as.numeric(table(factor(e_k, levels = E_star)))
  cat("  Method 1 (e_k only): sum(c_k) =", sum(c_k_method1), "\n")
  cat("    Expected: K_tilde =", data_list$cpp_data$K_tilde, "\n")
  
  # Method 2: count e_k and t_u (used in helper_functions.R)
  c_k_method2 <- as.numeric(table(factor(c(e_k, t_u), levels = E_star)))
  cat("  Method 2 (e_k + t_u): sum(c_k) =", sum(c_k_method2), "\n")
  cat("    Expected: K_tilde + U =", data_list$cpp_data$K_tilde + data_list$cpp_data$U, "\n")
  
  cat("  c_k used in data_list:", sum(data_list$cpp_data$c_k), "\n")
  cat("  Difference (method2 - method1):", sum(c_k_method2 - c_k_method1), "\n")
  
  # Set up initial values
  set.seed(1)
  z_init <- runif(data_list$cpp_data$I_mark)
  z_init <- z_init/sum(z_init)
  lambda_init <- rep(0.1, data_list$cpp_data$N)
  
  max_iter <- 1
  
  # Test C++ implementation
  cat("\n")
  cat("C++ Implementation:\n")
  res_cpp <- em_fit(make_model_data(data_list$cpp_data),
                    z_init = z_init,
                    lambda_init = lambda_init,
                    max_iter = max_iter, 
                    tol = 0, 
                    verbose = FALSE)
  
  cat("  sum(z_i) =", sum(res_cpp$z_i), "\n")
  cat("  sum(lambda_n) =", sum(res_cpp$lambda_n), "\n")
  
  # Test R implementation
  cat("\n")
  cat("R Implementation:\n")
  res_em <- do.call(
    em_estimate_raw,
    c(list(
      verbose = FALSE,
      max_iter = max_iter,
      tol = 0,
      z_init = z_init,
      lambda_init = lambda_init
    ), data_list$r_data)
  )
  
  cat("  sum(z) =", sum(res_em$z), "\n")
  cat("  sum(lambda) =", sum(res_em$lambda), "\n")
  
  # Detailed breakdown of z_i components
  cat("\n")
  cat("z_i Breakdown (R implementation):\n")
  cat("  z_head (i in 1:I): sum =", sum(res_em$z_head), "\n")
  cat("  z_tail (i in I+1:I_mark): sum =", sum(res_em$z_tail), "\n")
  
  # Check the components going into z calculation
  cat("\n")
  cat("Component Matrices:\n")
  cat("  sum(mu_mi) =", sum(res_em$mu_mi), "\n")
  cat("  sum(mu_bar_ji) =", sum(res_em$mu_bar_ji), "\n")
  cat("  sum(eta_ui) =", sum(res_em$eta_ui), "\n")
  cat("  sum(gamma_ci) =", sum(res_em$gamma_ci), "\n")
  
  cat("\n")
  cat("For z_head (equation 23):\n")
  numerator_head <- colSums(res_em$mu_mi) + 
                    colSums(res_em$mu_bar_ji[, 1:data_list$cpp_data$I]) + 
                    colSums(res_em$eta_ui[, 1:data_list$cpp_data$I]) + 
                    colSums(res_em$gamma_ci[, 1:data_list$cpp_data$I])
  cat("  sum(numerator) =", sum(numerator_head), "\n")
  cat("  N_star =", data_list$cpp_data$N_star, "\n")
  cat("  sum(z_head) = sum(numerator)/N_star =", sum(numerator_head)/data_list$cpp_data$N_star, "\n")
  
  if (data_list$cpp_data$K > 0) {
    cat("\n")
    cat("For z_tail (equation 24):\n")
    I <- data_list$cpp_data$I
    I_mark <- data_list$cpp_data$I_mark
    numerator_tail <- data_list$cpp_data$c_k + 
                      colSums(res_em$mu_bar_ji[, (I+1):I_mark]) + 
                      colSums(res_em$eta_ui[, (I+1):I_mark]) + 
                      colSums(res_em$gamma_ci[, (I+1):I_mark])
    cat("  sum(c_k) =", sum(data_list$cpp_data$c_k), "\n")
    cat("  sum(mu_bar for tail) =", sum(colSums(res_em$mu_bar_ji[, (I+1):I_mark])), "\n")
    cat("  sum(eta for tail) =", sum(colSums(res_em$eta_ui[, (I+1):I_mark])), "\n")
    cat("  sum(gamma for tail) =", sum(colSums(res_em$gamma_ci[, (I+1):I_mark])), "\n")
    cat("  sum(numerator) =", sum(numerator_tail), "\n")
    cat("  N_star =", data_list$cpp_data$N_star, "\n")
    cat("  sum(z_tail) = sum(numerator)/N_star =", sum(numerator_tail)/data_list$cpp_data$N_star, "\n")
  }
  
  cat("\n")
  cat("Total z sum check:\n")
  cat("  sum(z_head) + sum(z_tail) =", sum(res_em$z_head) + sum(res_em$z_tail), "\n")
  cat("  Expected: 1.0\n")
  
  # Return results for further inspection
  invisible(list(
    data_list = data_list,
    res_cpp = res_cpp,
    res_em = res_em,
    c_k_method1 = c_k_method1,
    c_k_method2 = c_k_method2
  ))
}

# ============================================================================
# Test Case 1: Small sample with constant hazards
# ============================================================================
cat("\n\n")
cat("############################################################\n")
cat("# TEST SUITE\n")
cat("############################################################\n")

set.seed(123)
sim_data_small <- simulate_idm_constant_hazards(
  n = 50, 
  a12 = 0.0008, 
  a13 = 0.0002, 
  a23 = 0.0016
)

result1 <- test_z_sum(sim_data_small, "Small sample (n=50) with constant hazards")

# ============================================================================
# Test Case 2: Medium sample with constant hazards
# ============================================================================
set.seed(456)
sim_data_medium <- simulate_idm_constant_hazards(
  n = 250, 
  a12 = 0.0008, 
  a13 = 0.0002, 
  a23 = 0.0016
)

result2 <- test_z_sum(sim_data_medium, "Medium sample (n=250) with constant hazards")

# ============================================================================
# Test Case 3: Larger sample
# ============================================================================
set.seed(789)
sim_data_large <- simulate_idm_constant_hazards(
  n = 500, 
  a12 = 0.001, 
  a13 = 0.0003, 
  a23 = 0.002
)

result3 <- test_z_sum(sim_data_large, "Large sample (n=500) with constant hazards")

# ============================================================================
# Test Case 4: Different hazard rates
# ============================================================================
set.seed(111)
sim_data_diff <- simulate_idm_constant_hazards(
  n = 200, 
  a12 = 0.002,  # Higher illness rate
  a13 = 0.0001, # Lower direct death
  a23 = 0.003   # Higher death after illness
)

result4 <- test_z_sum(sim_data_diff, "Different hazards (high a12, low a13)")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n\n")
cat("############################################################\n")
cat("# SUMMARY\n")
cat("############################################################\n")

cat("\nKey findings to check:\n")
cat("1. Does c_k differ between method 1 and method 2?\n")
cat("2. Is the sum(z_i) != 1 consistent across all test cases?\n")
cat("3. Does sum(numerator) in equations (23)+(24) equal N_star?\n")
cat("\n")
