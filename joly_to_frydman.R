df <- simulate_data(1000,censoring_time = 10)

data <- df$data$obs %>% arrange(status)
data

setup_frydman <- function(data) {
  
  data <- as.data.table(data)
  stopifnot(all(sort(names(data)) == sort(c("V_0", "V_healthy",
                                            "V_ill", "T_obs", "status"))))
  
  # case 1
  case_1 <- data[status == 1]
  case_1_exact <- case_1[V_healthy == T_obs]
  
  J <- case_1_exact[,.N]
  s_j <- case_1_exact$T_obs
  
  
  case_1_rest <- case_1[!(V_healthy == T_obs)]
  C <- case_1_rest[,.N]
  L_c <- case_1_rest$V_healthy
  t_c <- case_1_rest$T_obs
  
  # case 2
  case_2 <- data[status == 2]
  case_2_exact <- case_2[V_healthy == T_obs]
  K_tilde <- case_2_exact[,.N]
  e_k <- case_2_exact$T_obs
      
  case_2_rest <- case_2[!(V_healthy == T_obs)]
  U <- case_2_rest[,.N]  
  L_u <- case_2_rest$V_healthy
  t_u <- case_2_rest$T_obs
  
  # case 4
  case_4 <- data[status == 4]
  N_tilde <- case_4[,.N]  

  # case 3
  case_3_4 <- data[status %in% c(3,4)]
  M <- case_3_4[,.N]
  
  L_m <- case_3_4$V_healthy
  R_m <- case_3_4$V_ill
  t_m <- case_3_4$T_obs
  
  
  
  ##### K: E* - Obs and potential 1 -> 3 ####
  E_star <- unique(c(e_k, t_u))
  #c_k <- as.numeric(table(factor(c(e_k, t_u), levels = E_star)))
  # DANGER CHAT SAYS THAT THIS IS CORRECT
  # as sum(c_k) then is equal to K_tilde, otherwise K < sum(c_k) <= K + U, with = when e_k ∩ t_u = ∅
  c_k <- as.numeric(table(factor(e_k, levels = E_star)))
  K <- length(E_star)
  
  ##### N: T* - Obs and potential entry to state 3 from state 2: 1 -> 2 -> 3 ####
  T_star <- unique(c(t_m[in_N_tilde], t_u))
  d_n <- as.numeric(table(factor(c(t_m[in_N_tilde], t_u), levels = T_star)))
  
  N1_obs_of_T_star <- length(unique(t_m[in_N_tilde]))
  U_pos_obs_of_T_star <- length(setdiff(unique(t_u), unique(t_m[in_N_tilde])))
  
  N <- length(T_star)
  
  ##### Total: N* = M + U + C + K_tilde + J ####
  N_star <- M + U + C + K_tilde + J # Total count
  
  ##### Max 1 -> 2: M' = M + U + C ####
  M_mark <- M + U + C # Max number through 2.
  
  #### Creation of A sets ####
  
  ##### M: A_m := [L_m, R_m] ####
  A_m <- as.interval(matrix(c(L_m, R_m), ncol = 2, byrow = F))
  
  ##### W: M < m <= W := M + U, R_{M+u} = t_{M+u} ####
  A_u <- as.interval(matrix(c(L_u, t_u), ncol = 2, byrow = F))
  W = M + U
  
  ##### M': W := M + U < m <= M', R_{W+c} = t_{W+c} ####
  A_c <- as.interval(matrix(c(L_c, t_c), ncol = 2, byrow = F))
  
  ##### full_A_m: A_m ∪ A_u ∪ A_c ####
  full_A_m <- as.interval(rbind(A_m, A_u, A_c))
  
  ##### A := ⋃_{m=1}^{M'} A_m ####
  A_union <- get_interval(full_A_m)
  
  #### Data manipulation ####
  ##### I: Q_i = [l_i,r_i] ####
  
  # s_max = max(s_j, 1 <= j <= J)
  s_max <- max(s_j)
  
  # R_max = max(R_m, 1 <= m <= W)
  R_max <- max(A_m[, 2], A_u[, 2])
  
  # e*_max = max(e*_k, 1 <= k <= K)
  e_star_max <- max(E_star)
  
  # L_bar ={L_m, 1 <= m <= M'} ∪ {T* ∩ A} ∪ {S_J ∩ A} ∪ {s_max : s_max > R_max ∨ e*_max}
  L_bar <- c(
    full_A_m[, 1],
    intersect(A_union, T_star),
    intersect(A_union, s_j),
    na.omit(ifelse(s_max > max(R_max, e_star_max), s_max, NA))
  )
  
  # R_bar = {R_m, 1 <= m <= W} ∪ {∞}
  R_bar <- c(full_A_m[1:(M + U), 2], Inf)
  
  # !!!DANGER I AM UNSURE ABOUT THE CREATION OF Q!!!
  Q_i <- make_Q(L_bar, L_bar)
  
  Q <- get_interval(Q_i)
  
  I <- nrow(Q_i)
  
  
  ##### I'= K + I: Q_i' = e*_i-I ####
  Q_i_mark <- E_star
  
  ##### Q_full = list ####
  Q_full <- list(Q_i, Q_i_mark)
  
  ##### C: s_J+c = t_W+c ####
  s_j_c <- t_c
  
  s_j_full <- c(s_j, s_j_c)
  
  
  ##### N: lambda_n and I': z_i ####
  # Comment: I believe we have I' z_i's and N: lambda_n
  I_mark <- I + K
  
  
  list(
    Q_full   = Q_full,
    s_j_full = s_j_full, # for cal_alpha
    Q_i      = Q_i,
    full_A_m = full_A_m, # for cal_beta
    A_m      = A_m,
    A_u      = A_u,
    A_c      = A_c,      # for cal_* using M/U/C components
    T_star   = T_star,
    E_star   = E_star,   # for *star* arguments
    t_m      = t_m,
    t_u      = t_u,
    t_c      = t_c,      # event-time indices for M/U/C
    N_star   = N_star,   # denominator constant in (23)/(24)
    d_n      = d_n,      # first-term vector in (25); include zeros if not applicable
    c_k      = c_k       # constants for i > I in (24); optional
  )
}

list_1 <- setup_frydman(data)


res_em <- do.call(
  em_estimate,
  c(
    list(
      z_init = z_i,
      lambda_init = lambda_n,
      verbose = TRUE,
      max_iter = 200
    ),
    list1
  )
)
