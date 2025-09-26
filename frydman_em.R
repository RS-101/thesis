set.seed(1)
debugSource("frydman_setup.R")
debugSource("frydman_cal_est.R")
debugSource("frydman_plot.R")

# Replace NA with 0 (works for vectors and matrices)
na0 <- function(x) {
  x[is.na(x)] <- 0
  x
}

check_row_sums <- function(mu_mi, mu_bar_ji, eta_ui, gamma_ci, tol = 1e-6) {
  bad_mu    <- which(abs(rowSums(mu_mi)      - 1) >= tol)
  bad_mubar <- which(abs(rowSums(mu_bar_ji)  - 1) >= tol)
  bad_eta   <- which(abs(rowSums(eta_ui)     - 1) >= tol)
  bad_gamma <- which(abs(rowSums(gamma_ci)   - 1) >= tol)
  
  if (length(bad_mu))    cat("Rows of mu_mi not summing to 1:",    bad_mu,    "\n")
  if (length(bad_mubar)) cat("Rows of mu_bar_ji not summing to 1:", bad_mubar, "\n")
  if (length(bad_eta))   cat("Rows of eta_ui not summing to 1:",   bad_eta,   "\n")
  if (length(bad_gamma)) cat("Rows of gamma_ci not summing to 1:", bad_gamma, "\n")
  
  invisible(list(
    mu    = abs(rowSums(mu_mi)     - 1),
    mubar = abs(rowSums(mu_bar_ji) - 1),
    eta   = abs(rowSums(eta_ui)    - 1),
    gamma = abs(rowSums(gamma_ci)  - 1)
  ))
}


# EM to estimate (z, lambda) using Eqs. (23)–(25) and the provided helpers.
# Required inputs are exactly those needed by the helper functions.
em_estimate <- function(
    # Initial values
    z_init, lambda_init,
    # Data / design pieces used by the helper functions
    Q_full, s_j_full, # for cal_alpha
    Q_i, full_A_m, # for cal_beta
    A_m, A_u, A_c, # for cal_* using M/U/C components
    T_star, E_star, # for *star* arguments
    t_m, t_u, t_c, # event-time indices for M/U/C
    N_star, # denominator constant in (23)/(24)
    d_n, # first-term vector in (25); include zeros if not applicable
    c_k, # constants for i > I in (24); optional
    # Control
    max_iter = 200, tol = 1e-8, verbose = FALSE) {
  z_i <- as.numeric(z_init)
  lambda_n <- as.numeric(lambda_init)

  # Dimensions
  I <- length(z_i) # columns for the "I" block
  K <- if (is.null(c_k)) 0L else length(c_k) # optional tail (i > I)

  # History (optional)
  conv <- FALSE

  for (iter in seq_len(max_iter)) {
    z_prev <- z_i
    lambda_prev <- lambda_n

    ## ---------- E-step ----------
    alpha_ij <- cal_alpha(Q_full, s_j_full)
    beta_im <- cal_beta(Q_i, full_A_m)

    mu_mi <- cal_mu_MI(z_i, lambda_n, beta_im, Q_i, A_m, T_star)
    mu_mi <- na0(mu_mi)

    mu_bar_ji <- cal_mu_bar_JI_mark(alpha_ij, z_i)
    mu_bar_ji <- na0(mu_bar_ji)

    eta_ui <- cal_eta_UI_mark(z_i, lambda_n, beta_im, Q_i, A_u, E_star, T_star, t_u)
    eta_ui <- na0(eta_ui)

    gamma_ci <- cal_gamma_CI_mark(Q_i, A_c, t_c, T_star, lambda_n, alpha_ij, beta_im, z_i)
    gamma_ci <- na0(gamma_ci)

    rho_mn <- cal_rho_MN(t_m, T_star, A_m, mu_mi, Q_i)
    rho_mn <- na0(rho_mn)

    pi_un <- cal_pi_UN(t_u, T_star, A_u, eta_ui, Q_i)
    pi_un <- na0(pi_un)

    sigma_cn <- cal_sigma_CN(t_c, T_star, A_c, gamma_ci, Q_i)
    sigma_cn <- na0(sigma_cn)

    ## ---------- Verification -------
    # after computing mu_mi, mu_bar_ji, eta_ui, gamma_ci:
    check_row_sums(mu_mi, mu_bar_ji, eta_ui, gamma_ci)
    
    
    ## ---------- M-step ----------
    # z for 1..I via (23)
    z_head <- e_23(
      mu_mi          = mu_mi,
      mu_bar_ji      = mu_bar_ji,
      eta_ui         = eta_ui,
      gamma_ci       = gamma_ci,
      N_star         = N_star
    )

    # Optional tail i>I via (24)
    if (K > 0) {
      z_tail <- e_24(
        c_k          = c_k,
        mu_bar_ji    = mu_bar_ji,
        eta_ui       = eta_ui,
        gamma_ci     = gamma_ci,
        N_star       = N_star
      )
      z_i <- c(z_head, z_tail)
    } else {
      z_i <- z_head
    }

    # λ via (25)
    lambda_n <- e_25(
      d_n      = d_n,
      t_u      = t_u,
      T_star   = T_star,
      eta_ui   = eta_ui,
      rho_mn   = rho_mn,
      pi_un    = pi_un,
      sigma_cn = sigma_cn
    )


    ## ---------- Convergence ----------
    dz <- max(abs(z_i - z_prev))
    dl <- max(abs(lambda_n - lambda_prev))
    if (verbose) message(sprintf("iter %d: max|Δz|=%.3e, max|Δλ|=%.3e", iter, dz, dl))
    if (is.finite(dz) && is.finite(dl) && max(dz, dl) < tol) {
      conv <- TRUE
      break
    }
  }

  list(
    z = z_i,
    z_head = z_head,
    z_tail = z_tail,
    lambda = lambda_n,
    converged = conv,
    iterations = if (conv) iter else max_iter
  )
}

res_em <- em_estimate(z_init = z_i, lambda_init = lambda_n, Q_full, s_j_full, Q_i, full_A_m, A_m, A_u, A_c, T_star, E_star, t_m, t_u, t_c, N_star, d_n, c_k, verbose = T, max_iter = 200)

res_em

# We would expect sum z_p to be = 1
stopifnot(sum(res_em$z) - 1 < 1e-10)

estimators <- calc_F_and_hazards(
  grid_points = seq(0, 20, by = 0.5),
  z_i = res_em$z, lambda = res_em$lambda, Q_i, T_star, E_star
)

plot_estimators_gg(estimators)
