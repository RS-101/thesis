library(testthat)


setwd("/home/rasmus-emil/github/thesis/")
source("frydman/helper_functions_setup.R")
source("frydman/helper_function_generate_simple_data.R")
source("frydman/functions_em.R")
Rcpp::sourceCpp("frydman/functions_em.cpp")

test_that("product over t* matches between R and C++ implementations", {
  t_stars  <- seq(1, 2, by = 0.25)
  lambdas  <- (1:5) / 10
  
  R_res <- product_over_t_stars_one_interval(
    L = 0.9,
    R = 1.1,
    L_open = TRUE,
    R_open = TRUE,
    T_star = t_stars,
    lambda_n = lambdas
  )
  
  cpp_res <- product_over_t_stars_interval(
    L = 0.9,
    R = 1.1,
    L_open = TRUE,
    R_open = TRUE,
    T_star = t_stars,
    lambda_n = lambdas
  )
  
  # structure
  expect_equal(typeof(R_res), typeof(cpp_res))
  expect_equal(length(R_res), length(cpp_res))
  
  if (is.matrix(R_res) || is.array(R_res)) {
    expect_true(is.matrix(cpp_res) || is.array(cpp_res))
    expect_equal(dim(R_res), dim(cpp_res))
    expect_equal(rownames(R_res), rownames(cpp_res))
    expect_equal(colnames(R_res), colnames(cpp_res))
  } else {
    expect_false(is.matrix(cpp_res) || is.array(cpp_res))
  }
  
  # values
  expect_false(any(is.na(R_res)))
  expect_false(any(is.na(cpp_res)))
  expect_false(any(is.infinite(R_res)))
  expect_false(any(is.infinite(cpp_res)))
  expect_equal(c(R_res), c(cpp_res), tolerance = 1e-12)
  
  # order invariance (if functions are purely positional w.r.t. inputs)
  t_rev   <- rev(t_stars)
  lam_rev <- rev(lambdas)
  
  R_rev <- product_over_t_stars_one_interval(
    L = 0.9, R = 1.1, L_open = TRUE, R_open = TRUE,
    T_star = t_rev, lambda_n = lam_rev
  )
  cpp_rev <- product_over_t_stars_interval(
    L = 0.9, R = 1.1, L_open = TRUE, R_open = TRUE,
    T_star = t_rev, lambda_n = lam_rev
  )
  expect_equal(c(R_rev), c(cpp_rev), tolerance = 1e-12)
  
  # boundary variants (open/closed ends)
  for (L_open in c(TRUE, FALSE)) {
    for (R_open in c(TRUE, FALSE)) {
      R_b <- product_over_t_stars_one_interval(
        L = 0.9, R = 1.1, L_open = L_open, R_open = R_open,
        T_star = t_stars, lambda_n = lambdas
      )
      C_b <- product_over_t_stars_interval(
        L = 0.9, R = 1.1, L_open = L_open, R_open = R_open,
        T_star = t_stars, lambda_n = lambdas
      )
      expect_equal(c(R_b), c(C_b), tolerance = 1e-12)
    }
  }
})



test_that("R cal_* match C++ E-step pieces; (23)-(25) match M-step updates", {
  dat <- create_data(seed = 1)
  
  md_ptr <- setup_frydman_cpp(dat)

  I        <- nrow(dat$Q_i)
  I_mark   <- dat$I_mark
  N        <- length(dat$T_star)
  
  
  
  z_init      <- rep(1 / I_mark, I_mark)
  lambda_init <- rep(0.2, N)
  
  
  # ----- C++ run: one full EM iteration -----
  cpp <- em_fit(
    md_ptr,
    z_init = z_init,
    lambda_init = lambda_init,
    max_iter = 1,
    tol = 1e-12,
    verbose = FALSE
  )
  cpp_model_value <- model_data_to_list(md_ptr)
  
  # ----- R computations using the same inputs -----
  # alpha/beta
  alpha_R <- cal_alpha(dat$Q_i, dat$Q_i_mark, dat$s_j_full)
  beta_R  <- cal_beta(dat$Q_i, dat$full_A_m)
  
  # μ, mubar, η, γ (E-step with initial z/lambda)
  mu_R <- cal_mu_MI(z_init, lambda_init, beta_R, dat$Q_i, dat$A_m, dat$T_star);
  mubar_R <- cal_mu_bar_JI_mark(alpha_R, z_init, dat$J);
  eta_R <- cal_eta_UI_mark(z_init, lambda_init, beta_R, dat$Q_i, dat$A_u,
                           dat$E_star, dat$T_star, dat$t_u, dat$M)
  gamma_R <- cal_gamma_CI_mark(dat$Q_i, dat$A_c, dat$t_c, dat$T_star, lambda_init,
                               alpha_R, beta_R, z_init, dat$W, dat$J)
  
  # Aggregations for λ
  rho_R   <- cal_rho_MN(dat$t_m, dat$T_star, dat$A_m, mu_R, dat$Q_i)
  pi_R    <- cal_pi_UN(dat$t_u, dat$T_star, dat$A_u, eta_R, dat$Q_i)
  sigma_R <- cal_sigma_CN(dat$t_c, dat$T_star, dat$A_c, gamma_R, dat$Q_i)
  
  # M-step updates (23), (24), (25)
  z_head_R <- e_23(mu_R, mubar_R, eta_R, gamma_R, dat$N_star)
  if (length(dat$c_k) > 0) {
    z_tail_R <- e_24(dat$c_k, mubar_R, eta_R, gamma_R, dat$N_star, I = I, K_tilde = dat$K_tilde)
    z_R <- c(z_head_R, z_tail_R)
  } else {
    z_R <- z_head_R
  }
  lambda_R <- e_25(dat$d_n, dat$t_u, dat$T_star, eta_R, rho_R, pi_R, sigma_R, dat$N1_obs_of_T_star)
  # ----- Comparisons -----
  # alpha, beta
  expect_equal(dim(alpha_R), dim(cpp$alpha_ij))
  expect_equal(1*alpha_R, as.matrix(cpp$alpha_ij), tolerance = 1e-12, scale = 1)
  
  expect_equal(dim(beta_R), dim(cpp$beta_im))
  expect_equal(1*beta_R,  as.matrix(cpp$beta_im),  tolerance = 1e-12, scale = 1)
  # μ, mubar, η, γ
  expect_equal(dim(mu_R), dim(cpp$mu_mi))
  expect_equal(mu_R,    as.matrix(cpp$mu_mi),      tolerance = 1e-10, scale = 1)
  
  expect_equal(dim(mubar_R), dim(cpp$mu_bar_ji))
  expect_equal(mubar_R, as.matrix(cpp$mu_bar_ji),  tolerance = 1e-10, scale = 1)
  
  expect_equal(dim(eta_R), dim(cpp$eta_ui))
  expect_equal(eta_R,   as.matrix(cpp$eta_ui),     tolerance = 1e-10, scale = 1)
  
  expect_equal(dim(gamma_R), dim(cpp$gamma_ci))
  expect_equal(gamma_R, as.matrix(cpp$gamma_ci),   tolerance = 1e-10, scale = 1)
  
  # ρ, π, σ
  expect_equal(dim(rho_R), dim(cpp$rho_mn))
  expect_equal(rho_R,   as.matrix(cpp$rho_mn),     tolerance = 1e-12, scale = 1)
  
  expect_equal(dim(pi_R), dim(cpp$pi_un))
  expect_equal(pi_R,    as.matrix(cpp$pi_un),      tolerance = 1e-12, scale = 1)
  
  expect_equal(dim(sigma_R), dim(cpp$sigma_cn))
  expect_equal(sigma_R, as.matrix(cpp$sigma_cn),   tolerance = 1e-12, scale = 1)
  
  
  # z, λ after one M-step
  expect_equal(as.numeric(z_R),       as.numeric(cpp$z_i),       tolerance = 1e-10, scale = 1)
  expect_equal(as.numeric(lambda_R),  as.numeric(cpp$lambda_n),  tolerance = 1e-10, scale = 1)
  
  # Row-sum sanity checks (rows should sum to 1 within tolerance when defined)
  rs_mu    <- rowSums(mu_R)
  rs_mubar <- rowSums(mubar_R)
  rs_eta   <- rowSums(eta_R)
  if (nrow(gamma_R) > 0) {
    rs_gamma <- rowSums(gamma_R)
    expect_true(all(abs(rs_gamma - 1) < 1e-8 | (rs_gamma == 0)))
  }
  expect_true(all(abs(rs_mu - 1)    < 1e-8 | (rs_mu == 0)))
  expect_true(all(abs(rs_mubar - 1) < 1e-8 | (rs_mubar == 0)))
  expect_true(all(abs(rs_eta - 1)   < 1e-8 | (rs_eta == 0)))
})

