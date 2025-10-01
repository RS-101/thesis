# S3 method for plotting full_spline_hazard
# - cumulative=FALSE: plot hazards a12/a13/a23
# - cumulative=TRUE : plot cumulative hazards; for `x` use A12/A13/A23 if present,
#                     otherwise integrate its hazards; for `true_hazards` always integrate.
plot.full_spline_hazard <- function(x,
                                    true_hazards = NULL,
                                    xlim = c(0, 10),
                                    n = 256,
                                    cumulative = FALSE,
                                    title = if (cumulative) "Cumulative hazard" else "Hazard",
                                    ...) {
  stopifnot(is.numeric(xlim), length(xlim) == 2L, n >= 2)
  
  add_true <- !is.null(true_hazards)
  xvals <- seq(from = xlim[1], to = xlim[2], length.out = n)
  
  # Extract function lists -----------------------------------------------------
  as_fun_list <- function(obj) {
    # support either `obj$hazards` or the object itself being a list of functions
    f <- if (!is.null(obj$hazards)) obj$hazards else obj
    f[unlist(lapply(f, is.function))]
  }
  
  pick_named <- function(flist, pattern) {
    nm <- names(flist)
    if (is.null(nm)) stop("Function components must be named.")
    flist[grepl(pattern, nm)]
  }
  
  est_all <- as_fun_list(x)
  
  if (cumulative) {
    est_has_A <- any(grepl("^[A]", names(est_all)))
    est_fun <- if (est_has_A) pick_named(est_all, "^[A]") else pick_named(est_all, "^[a]")
    est_needs_integration <- !est_has_A
  } else {
    est_fun <- pick_named(est_all, "^[a]")
    est_needs_integration <- FALSE
  }
  
  if (add_true) {
    # true_hazards only contains hazards (lowercase)
    true_fun_haz <- pick_named(as_fun_list(true_hazards), "^[a]")
  }
  
  # Helpers to evaluate / integrate -------------------------------------------
  eval_on_grid <- function(f, grid) {
    vapply(grid, f, numeric(1))
  }
  
  integrate_on_grid <- function(f, grid, lower) {
    vapply(
      grid,
      function(t1) {
        if (t1 <= lower) return(0)
        stats::integrate(function(u) f(u), lower = lower, upper = t1, rel.tol = 1e-6)$value
      },
      numeric(1)
    )
  }
  
  build_df <- function(fun_list, type_label, cum = FALSE, lower = xlim[1]) {
    nms <- names(fun_list)
    rows <- lapply(seq_along(fun_list), function(i) {
      f <- fun_list[[i]]
      y <- if (cum) integrate_on_grid(f, xvals, lower) else eval_on_grid(f, xvals)
      # Use lowercase labels so A12 and a12 share colors across modes
      data.frame(x = xvals, y = y, fun = tolower(nms[i]), type = type_label,
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }
  
  # Data frames ----------------------------------------------------------------
  df_est <- if (cumulative && !est_needs_integration) {
    # Already cumulative functions (A12/A13/A23)
    build_df(est_fun, "estimate", cum = FALSE)
  } else if (cumulative && est_needs_integration) {
    # Integrate estimate hazards
    build_df(est_fun, "estimate", cum = TRUE, lower = xlim[1])
  } else {
    # Plain hazards
    build_df(est_fun, "estimate", cum = FALSE)
  }
  
  if (add_true) {
    df_true <- if (cumulative) {
      build_df(true_fun_haz, "true", cum = TRUE, lower = xlim[1])
    } else {
      build_df(true_fun_haz, "true", cum = FALSE)
    }
    df <- rbind(df_est, df_true)
  } else {
    df <- df_est
  }
  
  # Plot -----------------------------------------------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = x, y = y, color = fun, linetype = type)
  ) +
    ggplot2::geom_line(size = 1, alpha = 0.8) +
    ggplot2::scale_linetype_manual(values = c(estimate = "dashed", true = "solid")) +
    ggplot2::labs(
      title = title,
      x = "x", y = "value",
      color = "function", linetype = NULL
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::theme_minimal()
  
  if (!add_true) p <- p + ggplot2::guides(linetype = "none")
  p
}
