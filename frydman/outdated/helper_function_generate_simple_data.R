library(Rcpp)
source("frydman/helper_functions.R")
sourceCpp("optimizing_function_em.cpp")
#### Data generation #### 
create_data <- function(seed = NULL) { 
  if (!is.null(seed)) set.seed(seed)
  
  ## M : 1 -> 2 @ [L_m, R_m], potential 2 -> 3 at t_m (subset N_tilde observed)
  M <- 30
  N_tilde <- 15
  stopifnot(N_tilde <= M)
  
  id_m <- seq_len(M)
  in_N_tilde <- id_m %in% sample(id_m, size = N_tilde)
  
  L_m <- floor(runif(M, min = 1, max = 10))
  R_m <- L_m + runif(M, min = 1, max = 5)
  t_m <- R_m + runif(M, min = 0, max = 5)
  t_m_in_N_tilde <- t_m[in_N_tilde]
  
  ## K_tilde : 1 -> 3 @ e_k
  K_tilde <- 6
  e_k <- runif(K_tilde, min = 1, max = 10)
  
  ## J : 1 -> ? @ s_j
  J <- 30
  s_j <- runif(J, min = 1, max = 10)
  
  ## U : 1 @ L_u- -> ? -> 3 @ t_u
  U <- 24
  L_u <- runif(U, min = 1, max = 5)
  t_u <- L_u + runif(U, min = 2, max = 10)
  
  ## C : 1 @ L_c- -> ? -> ? @ t_c (censoring channel)
  C <- 10
  L_c <- runif(C, min = 1, max = 5)
  t_c <- L_c + runif(C, min = 2, max = 10)
  
  ## E* (observed/potential direct 1->3): union of e_k and t_u
  E_star <- sort(unique(c(e_k, t_u)))
  # c_k MUST sum to K_tilde; count e_k multiplicities only (not t_u)
  c_k <- as.numeric(table(factor(e_k, levels = E_star)))
  K <- length(E_star)
  
  ## T* (observed/potential entry to 3 from state 2): union of t_m[in_N_tilde] and t_u
  T_star <- sort(unique(c(t_m_in_N_tilde, t_u)))
  d_n <- as.numeric(table(factor(c(t_m_in_N_tilde, t_u), levels = T_star)))
  
  N1_obs_of_T_star <- length(unique(t_m_in_N_tilde))
  U_pos_obs_of_T_star <- length(setdiff(unique(t_u), unique(t_m_in_N_tilde)))
  N <- length(T_star)
  
  ## Totals
  N_star <- M + U + C + K_tilde + J
  M_mark <- M + U + C
  W <- M + U
  
  ## Build A-sets (2-col matrices). Keep them as plain matrices; helpers can wrap.
  A_m <- cbind(L_m, R_m)
  A_u <- cbind(L_u, t_u)
  A_c <- cbind(L_c, t_c)
  
  full_A_m <- rbind(A_m, A_u, A_c)
  
  ## A := ⋃_{m=1}^{M'} A_m, plus useful scalars
  # If helper_functions.R defines as.interval() and get_interval(), use them just for A_union.
  A_union <- get_interval(as.interval(full_A_m))
  
  s_max <- if (length(s_j)) max(s_j) else -Inf
  R_max <- max(c(A_m[, 2], A_u[, 2]), na.rm = TRUE)
  e_star_max <- if (length(E_star)) max(E_star) else -Inf
  
  ## Bars
  L_bar <- c(
    full_A_m[, 1],
    intersect(A_union, T_star),
    intersect(A_union, s_j),
    na.omit(ifelse(s_max > max(R_max, e_star_max), s_max, NA))
  )
  # Right endpoints for Q construction: all right bounds of [.,.) plus Inf
  R_bar <- c(full_A_m[seq_len(W), 2], Inf)
  
  ## Q_i : grid of closed intervals [l_i, r_i]
  # The previous code used make_Q(L_bar, L_bar) (bug). Use L_bar, R_bar.
  Q_i <- make_Q(L_bar, R_bar)
  stopifnot(is.matrix(Q_i), ncol(Q_i) == 2)
  
  Q <- get_interval(Q_i)
  I <- nrow(Q_i)
  
  ## Tail intervals / marks
  Q_i_mark <- E_star
  Q_full <- list(Q_i, Q_i_mark)
  
  s_j_c <- t_c
  s_j_full <- c(s_j, s_j_c)
  
  ## Initial values (not required by make_model_data, but useful for tests)
  I_mark <- I + K
  z_i <- runif(I_mark)
  lambda_n <- runif(N)
  
  model_data_list <- list(
    # ints
    J = J, C = C, K_tilde = K_tilde, U = U, N_tilde = N_tilde, M = M, W = W,
    N1_obs_of_T_star = N1_obs_of_T_star,
    U_pos_obs_of_T_star = U_pos_obs_of_T_star,
    N = N, N_star = N_star, M_mark = M_mark, I = I, K = K, I_mark = I_mark,
    
    # scalars
    s_max = s_max, R_max = R_max, e_star_max = e_star_max,
    
    # vectors
    s_j = s_j, L_c = L_c, t_c = t_c, e_k = e_k, L_u = L_u, t_u = t_u,
    t_m_in_N_tilde = t_m_in_N_tilde, L_m = L_m, R_m = R_m, t_m = t_m,
    E_star = E_star, T_star = T_star,
    L_bar = L_bar, R_bar = R_bar, s_j_c = s_j_c, s_j_full = s_j_full,
    Q = as.numeric(t(Q_i)), Q_i_mark = Q_i_mark, A_union = A_union,
    
    # integer vectors
    c_k = as.integer(c_k), d_n = as.integer(d_n),
    
    # matrices (2 columns)
    A_m = A_m, A_u = A_u, A_c = A_c,
    full_A_m = full_A_m, Q_i = Q_i
  )
  
  model_data_list_r <- list(
    Q_full = Q_full,
    s_j_full = s_j_full,
    Q_i = Q_i,
    full_A_m = full_A_m,
    A_m = A_m,
    A_u = A_u,
    A_c = A_c,
    T_star = T_star,
    E_star = E_star,
    t_m = t_m,
    t_u = t_u,
    t_c = t_c,
    N_star = N_star,
    d_n = d_n,
    c_k = c_k,
    I_mark = I_mark,
    J = J,
    M = M,
    W = W,
    K_tilde = K_tilde,
    N1_obs_of_T_star = N1_obs_of_T_star
  )
  

  list(cpp = model_data_list, 
       r = model_data_list_r)
}

my_data <- create_data()
md <- make_model_data(my_data$cpp)

cpp_fit <- em_fit(md)
# op <- options(warn = 2)      # treat warnings as errors
# try(em_fit(md))
# traceback()
# options(op)   

source("frydman/functions_em.R")

# Call with only valid args
r_fit <- do.call(
  em_estimate_raw,
  c(list(
    z_init      = rep(1/my_data$r$I_mark, my_data$r$I_mark),
    lambda_init = rep(1/2, length(my_data$r$T_star)),
    verbose     = TRUE,
    max_iter    = 1,
    tol         = 0.1
  ), my_data$r)
)

cpp_fit$z_i
r_fit$z





max(abs(as.numeric(cpp_fit$alpha_ij)-as.numeric(r_fit$alpha_ij)))
max(abs(as.numeric(cpp_fit$beta_im)-as.numeric(r_fit$beta_im)))
max(abs(as.numeric(cpp_fit$mu_mi)-as.numeric(r_fit$mu_mi)))
max(abs(as.numeric(cpp_fit$ws.mu_bar_ji)-as.numeric(r_fit$mu_bar_ji)))
max(abs(as.numeric(cpp_fit$ws.eta_ui)-as.numeric(r_fit$eta_ui)))
max(abs(as.numeric(cpp_fit$ws.gamma_ci)-as.numeric(r_fit$gamma_ci)))
max(abs(as.numeric(cpp_fit$z_i)-as.numeric(r_fit$z)))

