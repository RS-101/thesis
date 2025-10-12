source("joly/functions_likelihood.R")


# --- Packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(numDeriv)
  library(parallel)
})

# --- USER: supply these four functions -------------------------------------
# simulate_three_state(n, ...) -> data.frame with columns:
# V_0, V_healthy, V_ill, T_obs, status
# full_log_likehood(...) -> list(loglik = numeric(1))  # your existing function
# pen_mat_m_splines(knots) -> penalty matrix for M-splines
# make_all_3_spline_funs(theta, knots, degree) -> list of hazard functions
#
# Optionally, true hazards for accuracy scoring (ISE, etc.)
# h12_true(t), h13_true(t), h23_true(t)

# --- Fitter with approx-CV Vbar --------------------------------------------
do_likelihood_optim <- function(sim_dat, n_knots, degree = 3, penalizer, knots = NULL) {
  stopifnot(length(penalizer) == 3)
  n_par <- n_knots + degree - 1
  
  if (is.null(knots)) {
    knots_a12 <- seq(min(sim_dat$V_0), max(sim_dat$T_obs), length.out = n_knots)
    knots <- list(a12 = knots_a12, a13 = knots_a12, a23 = knots_a12)
  } else {
    stopifnot(length(knots$a12) == n_knots, length(knots$a13) == n_knots, length(knots$a23) == n_knots)
  }
  
  P12 <- pen_mat_m_splines(knots$a12)
  P13 <- pen_mat_m_splines(knots$a13)
  P23 <- pen_mat_m_splines(knots$a23)
  
  unpack_theta <- function(x) {
    list(
      a12 = x[1:n_par],
      a13 = x[(n_par+1):(2*n_par)],
      a23 = x[(2*n_par+1):(3*n_par)]
    )
  }
  
  ll_fun <- function(x) {
    theta <- unpack_theta(x)
    full_log_likehood(
      V_0       = sim_dat$V_0,
      V_healthy = sim_dat$V_healthy,
      V_ill     = sim_dat$V_ill,
      T_obs     = sim_dat$T_obs,
      status    = sim_dat$status,
      theta     = theta,
      degree    = degree,
      knots     = knots
    )$loglik
  }
  
  pen_fun <- function(x) {
    th <- unpack_theta(x)
    c(drop(t(th$a12) %*% P12 %*% th$a12),
      drop(t(th$a13) %*% P13 %*% th$a13),
      drop(t(th$a23) %*% P23 %*% th$a23))
  }
  
  pl_fun <- function(x) ll_fun(x) - sum(penalizer * pen_fun(x))
  
  obj_fun <- function(x) {
    val <- -pl_fun(x)
    if (!is.finite(val)) 1e10 else val
  }
  
  x0 <- rep(1, 3*(n_par))
  invisible(obj_fun(x0))
  
  res <- optim(
    par   = x0,
    fn    = obj_fun,
    method= "L-BFGS-B",
    lower = rep(1e-4, 3*n_par),
    upper = rep(13,   3*n_par)
  )
  
  theta_hat <- unpack_theta(res$par)
  pen_term_at_theta_hat <- sum(penalizer * pen_fun(res$par))
  est_hazards <- make_all_3_spline_funs(theta_hat, knots, degree)
  
  # Approximate CV score: Vbar = l(th_tilde) - tr(H_pen^{-1} H_ll)
  H_pen <- numDeriv::hessian(pl_fun, res$par)
  H_ll  <- numDeriv::hessian(ll_fun, res$par)
  
  # Stabilize inverse if needed
  make_pd <- function(M, eps = 1e-6, maxit = 5L) {
    i <- 0L
    while (i <= maxit) {
      ok <- try(chol(M), silent = TRUE)
      if (inherits(ok, "matrix")) break
      diag(M) <- diag(M) + eps
      eps <- eps * 10
      i <- i + 1L
    }
    M
  }
  H_pen <- make_pd(H_pen)
  tr_term <- sum(diag(solve(H_pen, H_ll)))
  Vbar    <- ll_fun(res$par) - tr_term
  
  list(
    loglik      = -res$value + pen_term_at_theta_hat,
    pen_loglik  = -res$value,
    approx_cv   = Vbar,
    trace_term  = tr_term,
    hazards     = est_hazards,
    theta       = theta_hat,
    knots       = knots,
    degree      = degree,
    optim       = res
  )
}

# --- Grid search over kappa (penalizer) ------------------------------------
fit_over_kappa <- function(sim_dat, n_knots, degree, kappa_grid, knots = NULL) {
  # kappa_grid is a data.frame with columns k12, k13, k23
  fits <- vector("list", nrow(kappa_grid))
  vbar <- numeric(nrow(kappa_grid))
  
  for (i in seq_len(nrow(kappa_grid))) {
    kap <- c(kappa_grid$k12[i], kappa_grid$k13[i], kappa_grid$k23[i])
    fit <- do_likelihood_optim(sim_dat, n_knots, degree, penalizer = kap, knots = knots)
    fits[[i]] <- fit
    vbar[i]   <- fit$approx_cv
  }
  best <- which.max(vbar)
  list(best_fit = fits[[best]],
       best_row = best,
       vbar     = vbar,
       fits     = fits)
}

# --- Accuracy metrics (example using ISE if true hazards are available) -----
ise_hazards <- function(h_est, h_true, t_grid) {
  mean((h_est(t_grid) - h_true(t_grid))^2) * (max(t_grid) - min(t_grid))
}

score_fit <- function(fit, t_grid, truths = NULL) {
  if (is.null(truths)) {
    return(list(approx_cv = fit$approx_cv))
  }
  h <- fit$hazards
  list(
    approx_cv = fit$approx_cv,
    ISE12 = ise_hazards(h$a12, truths$h12_true, t_grid),
    ISE13 = ise_hazards(h$a13, truths$h13_true, t_grid),
    ISE23 = ise_hazards(h$a23, truths$h23_true, t_grid)
  )
}

# --- One simulation replicate ----------------------------------------------
one_replicate <- function(seed,
                          n_train,
                          n_knots,
                          degree,
                          kappa_grid,
                          knots = NULL,
                          t_grid = NULL,
                          ...) {
  set.seed(seed)
  dat <- simulate_idm_weibull(1000,
                              shape12 = 1.1, scale12 = 1/0.0008,
                              shape13 = 1.8, scale13 = 1/0.0002,
                              shape23 = 1.3, scale23 = 1/0.0016)

  
  truths <- dat$true_data_generation$hazards
  dat <- dat$obs
  
  
  sel  <- fit_over_kappa(dat, n_knots, degree, kappa_grid, knots)
  fit  <- sel$best_fit
  out  <- score_fit(fit, t_grid = t_grid, truths = truths)
  c(list(seed = seed,
         best_row = sel$best_row,
         k12 = kappa_grid$k12[sel$best_row],
         k13 = kappa_grid$k13[sel$best_row],
         k23 = kappa_grid$k23[sel$best_row]),
    out)
}

# --- Simulation driver ------------------------------------------------------
run_simulation <- function(n_sims   = 200,
                           n_train  = 500,
                           n_knots  = 12,
                           degree   = 3,
                           k12_vals = 10^seq(-4, 2, length.out = 2),
                           k13_vals = 10^seq(-4, 2, length.out = 2),
                           k23_vals = 10^seq(-4, 2, length.out = 2),
                           knots    = NULL,
                           truths   = NULL,
                           t_grid   = seq(0, 5, length.out = 200),
                           ...      ) {
  
  kappa_grid <- expand.grid(k12 = k12_vals, k13 = k13_vals, k23 = k23_vals, KEEP.OUT.ATTRS = FALSE)
  seeds <- seq_len(n_sims)
  
  # parallel if available
  ncores <- max(1L, detectCores() - 1L)
  res_list <- mclapply(
    seeds,
    one_replicate,
    n_train    = n_train,
    n_knots    = n_knots,
    degree     = degree,
    kappa_grid = kappa_grid,
    knots      = knots,
    truths     = truths,
    t_grid     = t_grid,
    ...,
    mc.cores   = ncores
  )
  
  # bind results
  to_row <- function(x) {
    c(seed = x$seed,
      best_row = x$best_row,
      k12 = x$k12, k13 = x$k13, k23 = x$k23,
      approx_cv = x$approx_cv,
      ISE12 = if (!is.null(x$ISE12)) x$ISE12 else NA_real_,
      ISE13 = if (!is.null(x$ISE13)) x$ISE13 else NA_real_,
      ISE23 = if (!is.null(x$ISE23)) x$ISE23 else NA_real_)
  }
  M <- t(vapply(res_list, to_row, numeric(9)))
  out <- as.data.frame(M)
  
  summary <- within(list(), {
    kappa_selection_rate <- with(out, list(
      k12 = table(factor(k12, levels = k12_vals))/nrow(out),
      k13 = table(factor(k13, levels = k13_vals))/nrow(out),
      k23 = table(factor(k23, levels = k23_vals))/nrow(out)
    ))
    approx_cv_mean = mean(out$approx_cv, na.rm = TRUE)
    ISE_means = colMeans(out[, c("ISE12","ISE13","ISE23")], na.rm = TRUE)
  })
  
  list(results = out,
       grid    = kappa_grid,
       summary = summary)
}

# --- Example call (fill in your generator + truths) -------------------------
# truths <- list(h12_true = function(t) { ... },
#                 h13_true = function(t) { ... },
#                 h23_true = function(t) { ... })
sim_out <- run_simulation(
  n_sims  = 10,
  n_train = 40,
  n_knots = 10,
  t_grid  = seq(0, 6, length.out = 300)
)
str(sim_out$summary)
head(sim_out$results)
