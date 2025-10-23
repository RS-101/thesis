library(Rcpp)
library(data.table)
sourceCpp("optimizing_function_em.cpp")

source("datageneration/functions_simulate_data.R")
source("frydman/wrapper.R")
source("frydman/helper_functions.R")

# ============================================================================
# Quick validation check for z_i normalization
# Run this after any changes to the EM algorithm
# ============================================================================

quick_check <- function(n = 100, seed = 42) {
  set.seed(seed)
  
  cat("Running quick validation check (n =", n, ")...\n")
  
  # Simulate data
  sim_data <- simulate_idm_constant_hazards(
    n = n, 
    a12 = 0.0008, 
    a13 = 0.0002, 
    a23 = 0.0016
  )
  
  # Convert to list format
  data_list <- setup_data_to_list_format(sim_data$obs, add_r_format = FALSE)
  
  # Initialize
  set.seed(seed + 1)
  z_init <- runif(data_list$I_mark)
  z_init <- z_init/sum(z_init)
  lambda_init <- runif(data_list$N, min = 0.1, max = 0.5)
  
  # Run EM (single iteration)
  res_cpp <- em_fit(
    make_model_data(data_list),
    z_init = z_init,
    lambda_init = lambda_init,
    max_iter = 1, 
    tol = 0, 
    verbose = FALSE
  )
  
  # Check normalization
  z_sum <- sum(res_cpp$z_i)
  is_valid <- abs(z_sum - 1.0) < 1e-10
  
  # Report
  cat("\nResults:\n")
  cat("  Data: N_star =", data_list$N_star, 
      ", K_tilde =", data_list$K_tilde,
      ", U =", data_list$U, "\n")
  cat("  c_k: sum =", sum(data_list$c_k), 
      ", expected =", data_list$K_tilde,
      ", correct:", sum(data_list$c_k) == data_list$K_tilde, "\n")
  cat("  z_i: sum =", sprintf("%.15f", z_sum), "\n")
  cat("  Status:", ifelse(is_valid, "✓ PASS", "✗ FAIL"), "\n\n")
  
  if (!is_valid) {
    cat("WARNING: z_i does not sum to 1.0!\n")
    cat("Expected: 1.0, Got:", z_sum, "\n")
    cat("Difference:", z_sum - 1.0, "\n")
    stop("Validation check failed!")
  }
  
  invisible(list(z_sum = z_sum, valid = is_valid))
}

# Run the check
quick_check()

cat("========================================\n")
cat("Validation check PASSED! ✓\n")
cat("The z_i normalization is working correctly.\n")
cat("========================================\n")
