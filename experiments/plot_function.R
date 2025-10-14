library(splines2)

# Plots cumulative hazards A12, A13, A23.
# Color encodes transition (A12/A13/A23); linetype encodes source (True/NPMLE/PEN-MLE).
plot_cumhazards <- function(res, t_grid = NULL, t_max = NULL, n_grid = 400) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  library(ggplot2)
  
  # --- helpers ---------------------------------------------------------------
  .safe_get <- function(expr, default = NULL) {
    tryCatch(eval.parent(substitute(expr)), error = function(e) default)
  }
  .guess_tmax_from_stepfun <- function(f) {
    if (inherits(f, "stepfun")) {
      ks <- try(knots(f), silent = TRUE)
      if (!inherits(ks, "try-error") && length(ks)) return(max(ks, na.rm = TRUE))
    }
    NA_real_
  }
  .cum_from_hazard <- function(hfun, tvec) {
    if (!is.function(hfun)) return(rep(NA_real_, length(tvec)))
    vapply(tvec, function(tt) {
      if (tt <= 0) return(0)
      out <- try(integrate(function(u) hfun(u), lower = 0, upper = tt,
                           rel.tol = 1e-6, subdivisions = 100L), silent = TRUE)
      if (inherits(out, "try-error")) NA_real_ else out$value
    }, numeric(1))
  }
  .maybe_callable <- function(f, t_probe) {
    # Some objects are functions already, others are zero-arg factories returning a function
    if (is.function(f)) {
      test <- try(f(t_probe), silent = TRUE)
      if (!inherits(test, "try-error")) return(f)
      # maybe it is a factory needing '()'
      test2 <- try(f(), silent = TRUE)
      if (is.function(test2)) return(test2)
    }
    NULL
  }
  
  # --- pull true hazards (lowercase a**) -------------------------------------
  true_haz <- .safe_get(res$full_df$true_data_generation$hazards)
  a12 <- .safe_get(true_haz$a12)
  a13 <- .safe_get(true_haz$a13)
  a23 <- .safe_get(true_haz$a23)
  
  # --- pull estimators (cumulative hazards) ----------------------------------
  # NPMLE: typically Lambda**
  npm <- .safe_get(res$npmle$estimators, default = list())
  npm_A12 <- .safe_get(npm$Lambda12)
  npm_A13 <- .safe_get(npm$Lambda13)
  npm_A23 <- .safe_get(npm$Lambda23)
  
  # PEN-MLE: often in res$penmle$fit$hazards$A**
  pen_haz <- .safe_get(res$penmle$fit$hazards, default = list())
  pen_A12 <- .safe_get(pen_haz$A12)
  pen_A13 <- .safe_get(pen_haz$A13)
  pen_A23 <- .safe_get(pen_haz$A23)
  
  # --- choose time grid ------------------------------------------------------
  if (is.null(t_max)) {
    candidates <- c(
      .safe_get(max(res$full_df$obs$T_obs), NA_real_),
      .safe_get(res$full_df$true_data_generation$time_max, NA_real_),
      .safe_get(res$full_df$time_max, NA_real_),
      .safe_get(res$penmle$fit$time_max, NA_real_),
      .guess_tmax_from_stepfun(npm_A12),
      .guess_tmax_from_stepfun(npm_A13),
      .guess_tmax_from_stepfun(npm_A23)
    )
    t_max <- suppressWarnings(max(candidates, na.rm = TRUE))
    if (!is.finite(t_max) || is.na(t_max) || t_max <= 0) t_max <- 10
  }
  if (is.null(t_grid)) t_grid <- seq(0, t_max, length.out = n_grid)
  
  # Ensure estimator funcs are callable on vectors
  npm_A12 <- .maybe_callable(npm_A12, t_grid[1])
  npm_A13 <- .maybe_callable(npm_A13, t_grid[1])
  npm_A23 <- .maybe_callable(npm_A23, t_grid[1])
  pen_A12 <- .maybe_callable(pen_A12, t_grid[1])
  pen_A13 <- .maybe_callable(pen_A13, t_grid[1])
  pen_A23 <- .maybe_callable(pen_A23, t_grid[1])
  
  # --- build plotting data ---------------------------------------------------
  out_list <- list()
  
  # True (integrate lowercase hazards to get uppercase A**)
  out_list$true_A12 <- data.frame(
    t = t_grid, value = .cum_from_hazard(a12, t_grid),
    hazard = "A12", source = "True"
  )
  out_list$true_A13 <- data.frame(
    t = t_grid, value = .cum_from_hazard(a13, t_grid),
    hazard = "A13", source = "True"
  )
  out_list$true_A23 <- data.frame(
    t = t_grid, value = .cum_from_hazard(a23, t_grid),
    hazard = "A23", source = "True"
  )
  
  # NPMLE
  if (is.function(npm_A12)) out_list$npm_A12 <- data.frame(
    t = t_grid, value = as.numeric(suppressWarnings(npm_A12(t_grid))),
    hazard = "A12", source = "NPMLE"
  )
  if (is.function(npm_A13)) out_list$npm_A13 <- data.frame(
    t = t_grid, value = as.numeric(suppressWarnings(npm_A13(t_grid))),
    hazard = "A13", source = "NPMLE"
  )
  if (is.function(npm_A23)) out_list$npm_A23 <- data.frame(
    t = t_grid, value = as.numeric(suppressWarnings(npm_A23(t_grid))),
    hazard = "A23", source = "NPMLE"
  )
  
  # PEN-MLE
  if (is.function(pen_A12)) out_list$pen_A12 <- data.frame(
    t = t_grid, value = as.numeric(suppressWarnings(pen_A12(t_grid))),
    hazard = "A12", source = "PEN-MLE"
  )
  if (is.function(pen_A13)) out_list$pen_A13 <- data.frame(
    t = t_grid, value = as.numeric(suppressWarnings(pen_A13(t_grid))),
    hazard = "A13", source = "PEN-MLE"
  )
  if (is.function(pen_A23)) out_list$pen_A23 <- data.frame(
    t = t_grid, value = as.numeric(suppressWarnings(pen_A23(t_grid))),
    hazard = "A23", source = "PEN-MLE"
  )
  
  dat <- do.call(rbind, out_list)
  dat <- dat[is.finite(dat$value), ]
  
  p <- ggplot(dat, aes(x = t, y = value, color = hazard, linetype = source)) +
    geom_line(linewidth = 0.7) +
    labs(x = "Time", y = "Cumulative hazard",
         color = "Transition", linetype = "Type") +
    theme_minimal()
  
  print(p)
  invisible(list(data = dat, plot = p))
}
