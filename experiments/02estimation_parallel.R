# Parallelize with future.apply
library(qs)
library(future.apply)

# Use all data files (your original code only took the first)
files <- list.files("experiments/data/", full.names = TRUE)
dir.create("experiments/results", showWarnings = FALSE, recursive = TRUE)

# Optional: avoid CPU oversubscription if C++/BLAS code is multithreaded
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

# Plan workers (use all but one core)
future::plan(future::multisession, workers = max(1, future::availableCores() - 1))

run_one <- function(df_path) {
  # Load deps inside each worker
  library(qs)
  source("joly/functions_estimation.R")
  source("joly/functions_plot.R")
  source("frydman/wrapper.R")
  
  full_df  <- qread(df_path)
  obs_data <- full_df$obs
  
  est_npmle <- get_npmle(obs_data, max_iter = 100, tol = 0.001, verbose = TRUE)
  
  est_pl <- fit_idm(
    obs_data,
    n_knots = 7,
    degree = 3,
    kappa_values = c(0.1, 10, 100),
    verbose = TRUE
  )
  
  out_path <- file.path("experiments/results", paste0("estimate_", basename(df_path)))
  qsave(list(npmle = est_npmle, penmle = est_pl, full_df = full_df), file = out_path)
  out_path
}

# Parallel map; reproducible RNG across workers
res_paths <- future_lapply(files, run_one, future.seed = TRUE)
