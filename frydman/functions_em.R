source("frydman/helper_functions.R")

#### Helper function ####
##### Interval manipulation #####
###### intersect ######
# generic
intersect <- function(x, y, ...) {
  if (inherits(x, "interval")) {
    UseMethod("intersect", x)
  } else if (inherits(y, "interval")) {
    UseMethod("intersect", y)
  } else {
    base::intersect(x, y, ...)
  }
}

intersect.default <- function(x, y, ...) base::intersect(x, y, ...)

intersect.interval <- function(x, y) {
  if (inherits(y, "interval") & inherits(x, "numeric")) {
    x_temp <- x
    x <- y
    y <- x_temp
  }
  
  if(inherits(y, "numeric")) {
    L <- matrix(rep(x[,1], length(y)),ncol = nrow(x), byrow = T)
    R <- matrix(rep(x[,2], length(y)),ncol = nrow(x), byrow = T)
    
    if (inherits(x, "c_c")) {
      sel <- L <= y & y <= R
    } else if (inherits(x, "c_o")) {
      sel <- L <= y & y < R
    } else if (inherits(x, "o_c")) {
      sel <- L < y & y <= R
    } else if (inherits(x, "o_o")) {
      sel <- L < y & y < R
    }
    res = y[1 <= rowSums(sel)]
    return(res)
  } else if (inherits(y, "interval")) {
    stop("not implemented")
    if (inherits(x, "c_c")) {
      
    } else if (inherits(x, "c_o")) {
      
    } else if (inherits(x, "o_c")) {
      
    } else if (inherits(x, "o_o")) {
      
    }
  }
  stop("y should be numeric or interval")
}

###### get_interval ######

# generic
get_interval <- function(x, ...) {
  UseMethod("get_interval")
}

get_interval.matrix <- function(m, L_open = F, R_open = F) {
  if (ncol(m) != 2) stop("LR mat must be of dim m x 2")
  type = paste(ifelse(L_open, "o", "c"), ifelse(R_open, "o", "c"), sep = "_")
  if(nrow(m) > 1) {
    
    m <- unique(m)
    m_sorted <- m[order(m[, 1]), ]
    
    L <- m_sorted[,1][-1]
    R <- m_sorted[,2][-nrow(m)]
    
    
    if(L_open & R_open) {
      start_stop <- !(L < R)
    } else {
      start_stop <- !(L <= R)
    }
    interval_start <- c(min(m_sorted),L[start_stop])
    interval_end <- c(R[start_stop], max(m_sorted))
    res <- matrix(c(interval_start, interval_end),byrow = F, ncol = 2)
  } else {
    res <- m
  }
  class(res) <- c("interval", type, class(res))
  res
}


###### as.interval ######

as.interval <- function(x, L_open = F, R_open = F) {
  if (ncol(x) != 2) stop("LR mat must be of dim m x 2")
  type = paste(ifelse(L_open, "o", "c"), ifelse(R_open, "o", "c"), sep = "_")
  
  class(x) <- c("interval", type, class(x))
  x
}

###### contains ######

# generic
is_subset <- function(A, B, ...) {
  UseMethod("is_subset")
}

# checks if A ⊂ B
is_subset.interval <- function(A, B, strict_subset = F) {
  
  if(!(inherits(B, "interval") & all(dim(A) == dim(B)))) stop("B should be an of same dim as A interval")
  
  if (inherits(B, "c_c")) {
    l_compare <- `<=`
    r_compare <- `<=`
  } else if (inherits(B, "c_o")) {
    l_compare <- `<=`
    r_compare <- `<`
  } else if (inherits(B, "o_c")) {
    l_compare <- `<`
    r_compare <- `<=`
  } else if (inherits(B, "o_o")) {
    l_compare <- `<`
    r_compare <- `<`
  }
  
  if (inherits(A, "c_c")) {
  } else if (inherits(A, "c_o")) {
    r_compare <- `<=`
  } else if (inherits(A, "o_c")) {
    l_compare <- `<=`
  } else if (inherits(A, "o_o")) {
    l_compare <- `<=`
    r_compare <- `<=`
  }
  
  return(l_compare(B[,1],A[,1]) & r_compare(A[,2],B[,2]))
}




##### From paper specific #####

make_Q <- function(L_bar, R_bar) {
  L_bar <- sort(L_bar[!is.na(L_bar)])
  R_bar <- sort(R_bar[!is.na(R_bar)])
  Q <- matrix(c(rep(0L, length(L_bar)), rep(1L, length(R_bar)), 
                L_bar, R_bar), ncol = 2)
  Q <- Q[order(Q[,2]), ]
  tag <- which(diff(Q[, 1], 1) == 1)
  Q <- matrix(c(Q[tag, 2], Q[tag + 1, 2]), ncol = 2)
  
  Q <- as.interval(Q, L_open = F, R_open = F)
  Q
}

# Comment: the intervals input determines if the interval is open or closed
product_over_t_stars <- function(intervals, T_star, lambda_n) {
  if(!inherits(intervals, "interval")) stop("intervals need to be of type interval")
  T_stars_to_prod_over <- intersect(T_star, intervals)
  prod_lambdas <- lambda_n[which(T_star %in% T_stars_to_prod_over)]
  prod(1-prod_lambdas)
}

product_over_t_stars_one_interval <- function(L, R, L_open, R_open,T_star, lambda_n) {
  intervals <- as.interval(matrix(c(L, R), ncol = 2), L_open, R_open)
  T_stars_to_prod_over <- intersect(T_star, intervals)
  prod_lambdas <- lambda_n[which(T_star %in% T_stars_to_prod_over)]
  prod(1-prod_lambdas)
}






#### Prepare data frame to count data ####
prepare_data_R <- function(data) {
  
  data <- as.data.table(data)
  if(isFALSE("id" %in% names(data))) {data[, id:= .I]}
  
  stopifnot(all(sort(names(data)) == sort(c("id", "V_0", "V_healthy",
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
  t_m_in_N_tilde <- case_4$T_obs
  
  # case 3
  case_3 <- data[status == 3]
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
  T_star <- unique(c(t_m_in_N_tilde, t_u))
  d_n <- as.numeric(table(factor(c(t_m_in_N_tilde, t_u), levels = T_star)))
  
  N1_obs_of_T_star <- length(unique(t_m_in_N_tilde))
  U_pos_obs_of_T_star <- length(setdiff(unique(t_u), unique(t_m_in_N_tilde)))
  
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
    c_k      = c_k,       # constants for i > I in (24); optional,
    I_mark   = I_mark,
    J        = J, 
    M        = M,
    W        = W,
    K_tilde  = K_tilde,
    N1_obs_of_T_star = N1_obs_of_T_star
  )
}



#### Creating functions from (18) - (25) ####
# !!!DANGER SHOULD IT BE SHARP OR NOT!!!
# alpha is a I' x (J+C) matrix ~ alpha_ij = I(Q_i_full subset [s_j_full,inf)).

cal_alpha <- function(Q_i, Q_i_mark, s_j_full) {
  I <- nrow(Q_i)
  I_mark <- I + length(Q_i_mark)
  J_p_C <- length(s_j_full)
  alpha <- matrix(, nrow = I_mark, ncol = J_p_C)
  for (i in 1:I_mark) {
    if (i <= I) {
      alpha[i, ] <- s_j_full <= Q_i[i, 1]
    } else {
      alpha[i, ] <- s_j_full <= Q_i_mark[i - I]
    }
  }
  alpha
}

# !!!DANGER SHOULD IT BE SHARP OR NOT!!!
# beta is a I x M' matrix ~ beta_im = I(Q_i subset A_m).
cal_beta <- function(Q_i, full_A_m) {
  I <- nrow(Q_i)
  M_mark <- nrow(full_A_m)
  
  beta <- matrix(numeric(I * M_mark), nrow = I, ncol = M_mark)
  for (i in 1:I) {
    beta[i, ] <- full_A_m[, 1] <= Q_i[i, 1] & Q_i[i, 2] <= full_A_m[, 2]
  }
  beta
}


##### μ_mi(z,λ,β,Q,A_m) ∈ M x I, ####
cal_mu_MI <- function(z_i, lambda_n, beta_im, Q_i, A_m, T_star) {
  I <- nrow(Q_i)
  M <- nrow(A_m)
  
  beta_im <- beta_im[1:I, 1:M]
  if (!all(dim(beta_im) == c(I, M))) stop("beta wrong dim")
  
  mu <- matrix(nrow = M, ncol = I)
  
  r_i <- Q_i[, 2]
  R <- A_m[, 2]
  
  for (m in 1:M) {
    prod_res <- unlist(lapply(r_i, function(r_i) {
      product_over_t_stars_one_interval(
        r_i, R[m],
        L_open = TRUE,
        R_open = FALSE,
        T_star = T_star,
        lambda_n = lambda_n
      )
    }))
    denum <- sum(beta_im[, m] * z_i[1:I] * prod_res)
    for (i in 1:I) {
      num <- beta_im[i, m] * z_i[i] * prod_res[i]
      mu[m, i] <- ifelse(num != 0, num / denum, 0)
    }
  }
  mu
}


##### μ_bar_mi(z, ɑ) ∈ J x I', ####
cal_mu_bar_JI_mark <- function(alpha_ij, z_i, J) {
  I_mark <- nrow(alpha_ij)
  
  prod <- alpha_ij[1:I_mark, 1:J] * z_i
  res <- t(sweep(prod, 2, colSums(prod), "/"))
  res[is.na(res)] <- 0
  res
}

##### η_ji(z,λ,β,Q,A_m) ∈ U x I', ####
cal_eta_UI_mark <- function(z_i, lambda_n, beta_im, Q_i, A_u, E_star, T_star, t_u, M) {
  # Map between t_u and λ_n
  U <- nrow(A_u)
  I <- nrow(beta_im)
  I_mark <- length(z_i)
  
  eta_ui <- matrix(nrow = U, ncol = I_mark)
  
  r_i <- Q_i[, 2]
  for (u in 1:U) {
    lambda_M_p_u <- lambda_n[T_star == t_u[u]]
    prod_res <- unlist(lapply(r_i, function(r_i) {
      product_over_t_stars_one_interval(
        r_i, t_u[u],
        L_open = TRUE,
        R_open = TRUE,
        T_star = T_star,
        lambda_n = lambda_n
      )
    }))
    
    denom <- lambda_M_p_u * sum(prod_res * beta_im[1:I, (M + u)] * z_i[1:I]) +
      sum((E_star == t_u[u]) * z_i[(I + 1):I_mark])
    
    for (i in 1:I_mark) {
      if (i <= I) {
        eta_ui[u, i] <- lambda_M_p_u *
          prod_res[i] *
          beta_im[i, (M + u)] *
          z_i[i] / denom
      } else if (i <= I_mark) {
        eta_ui[u, i] <- as.numeric(t_u[u] == E_star[i - I]) * z_i[i] / denom
      }
    }
  }
  eta_ui
}

# eta_ui <- cal_eta_UI_mark(z_i, lambda_n, beta_im, Q_i, A_u, E_star, T_star, t_u)

##### γ_ji() ∈ C x I', ####
cal_gamma_CI_mark <- function(Q_i, A_c, t_c, T_star, lambda_n, alpha_ij, beta_im, z_i, W, J) {
  I <- nrow(Q_i)
  C <- nrow(A_c)
  I_mark <- length(z_i)
  gamma_ci <- matrix(nrow = C, ncol = I_mark)
  r_i <- Q_i[, 2]
  
  if (C == 0) return(gamma_ci)
  for (c in 1:C) {
    prod_res <- unlist(lapply(r_i, function(r_i) {
      product_over_t_stars_one_interval(
        r_i, t_c[c],
        L_open = TRUE,
        R_open = FALSE,
        T_star = T_star,
        lambda_n = lambda_n
      )
    }))
    
    denum <- sum(prod_res * beta_im[1:I, (W + c)] * z_i[1:I]) + sum(alpha_ij[, (J + c)] * z_i)      
    for (i in 1:I_mark) {
      if (i <= I) {
        gamma_ci[c,i] <- alpha_ij[i,(J+c)]*z_i[i]/denum +
          (prod_res[i]*beta_im[i,(W+c)]*z_i[i])/denum
      } else if (i <= I_mark) {
        gamma_ci[c,i] <- alpha_ij[i,(J+c)]*z_i[i]/denum
      }
    }
  }
  
  gamma_ci
}


##### ρ ∈ M x N ####
cal_rho_MN <- function(t_m, T_star, A_m, mu_mi, Q_i) {
  M <- length(t_m)
  N <- length(T_star)
  I <- nrow(Q_i)
  rho_mn <- matrix(nrow = M, ncol = N)
  
  for (n in 1:N) {
    for (m in 1:M) {
      if (t_m[m] >= T_star[n]) {
        int_L_t <- as.interval(
          matrix(rep(c(A_m[m,1], T_star[n]), I), ncol = 2, byrow = TRUE),
          L_open = FALSE, R_open = TRUE
        )
        rho_mn[m,n] <- sum(mu_mi[m,1:I] * is_subset(Q_i, int_L_t))
      } else {
        rho_mn[m,n] <- 0
      }
    }
  }  
  rho_mn
}

##### π ∈ U x N ####
cal_pi_UN <- function(t_u, T_star, A_u, eta_ui, Q_i) {
  U <- length(t_u)
  N <- length(T_star)
  I <- nrow(Q_i)
  pi_un <- matrix(nrow = U, ncol = N)
  
  for (n in 1:N) {
    for (u in 1:U) {
      if (t_u[u] >= T_star[n]) {
        int_L_t <- as.interval(
          matrix(rep(c(A_u[u,1], T_star[n]), I), ncol = 2, byrow = TRUE),
          L_open = FALSE, R_open = TRUE
        )
        pi_un[u,n] <- sum(eta_ui[u,1:I] * is_subset(Q_i, int_L_t))
      } else {
        pi_un[u,n] <- 0
      }
    }
  }  
  pi_un
}


# pi_un <- cal_pi_UN(t_u, T_star, A_u, eta_ui, Q_i)
##### σ ∈ C x N ####
cal_sigma_CN <- function(t_c, T_star, A_c, gamma_ci, Q_i) {
  C <- length(t_c)
  N <- length(T_star)
  I <- nrow(Q_i)
  
  sigma_cn <- matrix(nrow = C, ncol = N)
  
  for (n in 1:N) {
    for (c in 1:C) {
      if (t_c[c] >= T_star[n]) {
        int_L_t <- as.interval(
          matrix(rep(c(A_c[c,1], T_star[n]), I), ncol = 2, byrow = TRUE),
          L_open = FALSE, R_open = TRUE
        )
        sigma_cn[c,n] <- sum(gamma_ci[c,1:I] * is_subset(Q_i, int_L_t))
      } else {
        sigma_cn[c,n] <- 0
      }
    }
  }  
  sigma_cn
}
# sigma_cn <- cal_sigma_CN(t_c, T_star, A_c, gamma_ci,Q_i)

#### (23)-(25) ####
# (23)  z_i = (Σ_m μ_{mi} + Σ_j \bar{μ}_{ji} + Σ_u η_{ui} + Σ_c γ_{ci}) / N*
e_23 <- function(mu_mi, mu_bar_ji, eta_ui, gamma_ci, N_star) {
  stopifnot(is.matrix(mu_mi), is.matrix(mu_bar_ji),
            is.matrix(eta_ui), is.matrix(gamma_ci),
            length(N_star) == 1, is.finite(N_star))
  
  I <- ncol(mu_mi)
  mu_mi <- mu_mi[,1:I]
  mu_bar_ji <- mu_bar_ji[,1:I]
  eta_ui <- eta_ui[,1:I]
  gamma_ci <- gamma_ci[,1:I]
  mu_mi <- mu_mi[,1:I]
  
  if (!all(ncol(mu_mi) == I, ncol(mu_bar_ji) == I, ncol(eta_ui) == I, ncol(gamma_ci) == I)) {
    stop("All matrices must have the same number of columns (I).")
  }
  
  (colSums(mu_mi) + colSums(mu_bar_ji) + colSums(eta_ui) + colSums(gamma_ci)) / N_star
}

# (24)  for i > I: z_i = (c_{i-I} + Σ_j \bar{μ}_{ji} + Σ_u η_{ui} + Σ_c γ_{ci}) / N*
# Here c_vec is the vector (c_{1},…,c_{K}) for those extra indices,
# and the matrices contain only the corresponding K columns.
e_24 <- function(c_k, mu_bar_ji, eta_ui, gamma_ci, N_star, I, K_tilde) {
  stopifnot(is.numeric(c_k), is.matrix(mu_bar_ji),
            is.matrix(eta_ui), is.matrix(gamma_ci),
            length(N_star) == 1, is.finite(N_star))
  
  I_mark <- ncol(mu_bar_ji)
  I_diff <- I_mark - I
  mu_bar_ji <- mu_bar_ji[,(I+1):I_mark]
  eta_ui <- eta_ui[,(I+1):I_mark]
  gamma_ci <- gamma_ci[,(I+1):I_mark]
  
  if (!all(sum(c_k) == K_tilde, ncol(mu_bar_ji) == I_diff, ncol(eta_ui) == I_diff, ncol(gamma_ci) == I_diff)) {
    stop("mu_bar_ji, eta_ui, gamma_ci must each have K columns matching length(c_vec).")
  }
  
  (c_k + colSums(mu_bar_ji) + colSums(eta_ui) + colSums(gamma_ci)) / N_star
}

# λ_n from (25):
# λ_n = [ I(n ≤ N1) d_n + Σ_u I(t_{M+u} = t*_n) Σ_i η_{ui} ] /
#       [ Σ_m ρ_{mn} + Σ_u π_{un} + Σ_c σ_{cn} ]

e_25 <- function(d_n, t_u, T_star, eta_ui, rho_mn, pi_un, sigma_cn, N1_obs_of_T_star) {
  N <- length(T_star)
  stopifnot(is.numeric(d_n), length(d_n) > 0)
  #N <- length(d_n)
  stopifnot(length(T_star) == N,
            ncol(rho_mn) == N,
            ncol(pi_un) == N,
            nrow(eta_ui) == length(t_u))
  
  # Denominator: column sums across m, u, c
  denom <- colSums(rho_mn) + colSums(pi_un) + colSums(sigma_cn)
  
  # Numerator second term: for each n, sum_u 1(t_{M+u} = t*_n) * Σ_i η_{ui}
  s_eta_u   <- rowSums(eta_ui)                              # Σ_i η_{ui}
  s_by_time <- rowsum(s_eta_u, group = as.character(t_u))
  s_by_time <- setNames(as.numeric(s_by_time), rownames(s_by_time))
  eta_term  <- s_by_time[as.character(T_star)]
  eta_term[is.na(eta_term)] <- 0
  
  
  # Numerator first term: I(n ≤ N1) d_n
  num <- as.numeric(seq_len(N) <= N1_obs_of_T_star) * d_n + eta_term
  
  num / denom
}

set.seed(1)

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
#### Run EM function ####

# EM to estimate (z, lambda) using Eqs. (23)–(25) and the provided helpers.
# Required inputs are exactly those needed by the helper functions.
em_estimate_raw <- function(
    # Initial values
  z_init = NULL, lambda_init = NULL,
  # Data / design pieces used by the helper functions
  Q_full, s_j_full, # for cal_alpha
  Q_i, full_A_m, # for cal_beta
  A_m, A_u, A_c, # for cal_* using M/U/C components
  T_star, E_star, # for *star* arguments
  t_m, t_u, t_c, # event-time indices for M/U/C
  N_star, # denominator constant in (23)/(24)
  d_n, # first-term vector in (25); include zeros if not applicable
  c_k, # constants for i > I in (24); optional
  I_mark,
  J,
  M,
  W,
  K_tilde,
  N1_obs_of_T_star,
  # Control
  max_iter = 200, tol = 1e-8, verbose = FALSE) {
  
  I <- nrow(Q_i)
  N <- length(T_star)
  
  if(is.null(z_init) | is.null(lambda_init)) {
    z_init <- runif(I_mark)
    z_init <- z_init/sum(z_init)
    
    lambda_init <- runif(N)
  }
  z_i <- as.numeric(z_init)
  lambda_n <- as.numeric(lambda_init)
  
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
    
    mu_bar_ji <- cal_mu_bar_JI_mark(alpha_ij, z_i, J)
    mu_bar_ji <- na0(mu_bar_ji)
    
    eta_ui <- cal_eta_UI_mark(z_i, lambda_n, beta_im, Q_i, A_u, E_star, T_star, t_u, M)
    eta_ui <- na0(eta_ui)
    
    gamma_ci <- cal_gamma_CI_mark(Q_i, A_c, t_c, T_star, lambda_n, alpha_ij, beta_im, z_i, W)
    gamma_ci <- na0(gamma_ci)
    
    rho_mn <- cal_rho_MN(t_m, T_star, A_m, mu_mi, Q_i)
    rho_mn <- na0(rho_mn)
    
    pi_un <- cal_pi_UN(t_u, T_star, A_u, eta_ui, Q_i)
    pi_un <- na0(pi_un)
    
    if(length(t_c) > 0){
      sigma_cn <- cal_sigma_CN(t_c, T_star, A_c, gamma_ci, Q_i)
      sigma_cn <- na0(sigma_cn)
    } else {
      sigma_cn = as.matrix(0)
    }
    
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
        N_star       = N_star,
        I            = I,
        K_tilde      = K_tilde
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
      sigma_cn = sigma_cn,
      N1_obs_of_T_star = N1_obs_of_T_star
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



#### Wrapper ####
get_npmle_r <- function(data, max_iter = 200, tol = 1e-8, verbose = FALSE) {
  list_data <- prepare_data_R(data = data)
  
  res_em <- do.call(
    em_estimate_raw,
    c(list(
      verbose = verbose,
      max_iter = max_iter,
      tol = tol
    ),list_data)
  )
  
  estimators <- calc_F_and_hazards(
    grid_points = seq(0, 3, by = 0.01),
    z_i = res_em$z, lambda = res_em$lambda, 
    list_1$Q_full, list_1$T_star, list_1$E_star
  )
  
  plot <- plot_estimators_gg(estimators)
  
  
  list(
    estimators = estimators,
    plot = plot, 
    settings = list(
      data = data,
      verbose = verbose,
      max_iter = max_iter,
      tol = tol
    ))
}
  





