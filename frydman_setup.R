
rm(list = ls())
debugSource("aux_frydman.R")

#### Data generation ####
##### M : 1 -> 2 @ L,R -> 3 or N_tilde : ? @ t_m  ####
M <- 30 # interval transitions 1 -> 2.
N_tilde <- 15 # subset of 1:M 1 -> 2 -> 3 the rest of M is censored at 2.
if (N_tilde > M) stop("N_tilde should be lower than M")

id_m <- seq(M)
in_N_tilde <- id_m %in% sample(id_m, size = N_tilde)

L_m <- floor(runif(M, min = 1, max = 10))
R_m <- L_m + floor(runif(M, min = 1, max = 5))

t_m <- R_m + floor(runif(M, min = 0, max = 5))


##### K_tilde : 1 -> 3 @ e_k ####
K_tilde <- 6 # 1 -> 3
e_k <- floor(runif(K_tilde, min = 1, max = 10))

##### J : 1 -> ? @ s_j ####
J <- 30
s_j <- floor(runif(J, min = 1, max = 10))

##### U: 1 @ L_u- -> ? -> 3 @ t_u  ####
U <- 24
L_u <- floor(runif(U, min = 1, max = 5)) # last time in state 1
t_u <- L_u + floor(runif(U, min = 2, max = 10))

##### C: 1 @ L_c- -> ? -> ? @ t_c  ####
C <- 10
L_c <- floor(runif(C, min = 1, max = 5)) # last time in state 1
t_c <- L_c + floor(runif(C, min = 2, max = 10)) # censored

##### K: E* - Obs and potential 1 -> 3 ####
E_star <- unique(c(e_k, t_u))
c_k <- as.numeric(table(factor(c(e_k, t_u), levels = E_star)))
K <- length(E_star)

##### N: T* - Obs and potential 1 -> 2 -> 3 ####
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
A_union <- get_interval(as.interval(matrix(c(1, 4), nrow = 1, byrow = T)))

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
  intersect(s_j, A_union),
  na.omit(ifelse(s_max > max(R_max, e_star_max), s_max, NA))
)

# R_bar = {R_m, 1 <= m <= W} ∪ {∞}
R_bar <- c(full_A_m[1:(M + U), 2], Inf)

L_bar_sort_distinct <- sort(unique(L_bar))
R_bar_sort_distinct <- sort(unique(R_bar))

# !!!DANGER I AM UNSURE ABOUT THE CREATION OF Q!!!
Q_i <- make_Q(L_bar_sort_distinct, R_bar_sort_distinct)

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

z_i <- runif(I_mark)
lambda_n <- runif(N)

#### Creating functions from (18) - (25) ####
# !!!DANGER SHOULD IT BE SHARP OR NOT!!!
# alpha is a I' x (J+C) matrix ~ alpha_ij = I(Q_i_full subset [s_j_full,inf)).

cal_alpha <- function(Q_full, s_j_full) {
  Q_i <- Q_full[[1]]
  Q_i_mark <- Q_full[[2]]
  alpha <- matrix(, nrow = I_mark, ncol = J + C)
  for (i in 1:I_mark) {
    if (i <= I) {
      alpha[i, ] <- s_j_full <= Q_i[i, 1]
    } else {
      alpha[i, ] <- s_j_full <= Q_i_mark[i - I]
    }
  }
  alpha
}

alpha_ij <- cal_alpha(Q_full, s_j_full)
# !!!DANGER SHOULD IT BE SHARP OR NOT!!!
# beta is a I x M' matrix ~ beta_im = I(Q_i subset A_m).
cal_beta <- function(Q_i, full_A_m) {
  beta <- matrix(numeric(I * M_mark), nrow = I, ncol = M_mark)
  for (i in 1:I) {
    beta[i, ] <- full_A_m[, 1] <= Q_i[i, 1] & Q_i[i, 2] <= full_A_m[, 2]
  }
  beta
}
beta_im <- cal_beta(Q_i, full_A_m)


##### μ_mi(z,λ,β,Q,A_m) ∈ M x I, ####
cal_mu_MI <- function(z_i, lambda_n, beta_im, Q_i, A_m, T_star) {
  beta_im <- beta_im[1:I, 1:M]
  if (!all(dim(beta_im) == c(I, M))) stop("beta wrong dim")

  mu <- matrix(nrow = M, ncol = I)

  r_i <- Q_i[, 2]
  R <- A_m[, 2]

  for (m in 1:M) {
    prod_res <- unlist(lapply(r_i, function(r_i) {
      product_over_t_stars(
        intervals = as.interval(
          matrix(c(r_i, R[m]), ncol = 2),
          L_open = T,
          R_open = F
        ),
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
# mu_mi <- cal_mu_MI(z_i, lambda_n, beta_im, Q_i, A_m, T_star)

##### μ_bar_mi(z, ɑ) ∈ J x I', ####
cal_mu_bar_JI_mark <- function(alpha_ij, z_i) {
  prod <- alpha_ij * z_i
  res <- t(sweep(prod, 2, colSums(prod), "/"))
  res[is.na(res)] <- 0
  res
}

# mu_bar_ji <- cal_mu_bar_JI_mark(alpha_ij, z_i)

##### η_ji(z,λ,β,Q,A_m) ∈ U x I', ####
cal_eta_UI_mark <- function(z_i, lambda_n, beta_im, Q_i, A_u, E_star, T_star, t_u) {
  # Map between t_u and λ_n

  eta_ui <- matrix(nrow = U, ncol = I_mark)

  r_i <- Q_i[, 2]
  L_u_U <- A_u[, 1]^U
  for (u in 1:U) {
    lambda_M_p_u <- lambda_n[T_star == t_u[u]]
    prod_res <- unlist(lapply(r_i, function(r_i) {
      product_over_t_stars(
        intervals = as.interval(
          matrix(c(r_i, t_u[u]), ncol = 2),
          L_open = T,
          R_open = T
        ),
        T_star = T_star,
        lambda_n = lambda_n
      )
    }))
    for (i in 1:I_mark) {
      if (i <= I) {
        eta_ui[u, i] <- lambda_M_p_u *
          prod_res[i] *
          beta_im[i, (M + u)] *
          z_i[i] / L_u_U[u]
      } else if (i <= I_mark) {
        eta_ui[u, i] <- as.numeric(t_u[u] == E_star[i - I]) * z_i[i] / L_u_U[u]
      }
    }
  }
  eta_ui
}

# eta_ui <- cal_eta_UI_mark(z_i, lambda_n, beta_im, Q_i, A_u, E_star, T_star, t_u)

##### γ_ji() ∈ C x I', ####
cal_gamma_CI_mark <- function(Q_i, A_c, t_c, T_star, lambda_n, alpha_ij, beta_im, z_i) {
  gamma_ci <- matrix(nrow = C, ncol = I_mark)

  r_i <- Q_i[, 2]
  L_c_C <- A_c[, 1]^C
  
  for (c in 1:C) {
    prod_res <- unlist(lapply(r_i, function(r_i) {
      product_over_t_stars(
        intervals = as.interval(
          matrix(c(r_i, t_c[c]), ncol = 2),
          L_open = T,
          R_open = F
        ),
        T_star = T_star,
        lambda_n = lambda_n
      )
    }))
    for (i in 1:I_mark) {
      if (i <= I) {
        gamma_ci[c,i] <- alpha_ij[i,(J+c)]/L_c_C[c]+
          (prod_res[i]*beta_im[i,(W+c)])/L_c_C[c]
      } else if (i <= I_mark) {
        gamma_ci[c,i] <- alpha_ij[i,(J+c)]*z_i[i]/L_c_C[c]
      }
    }
  }

  gamma_ci
}

# gamma_ci <- cal_gamma_CI_mark(Q_i, A_c, t_c, T_star, lambda_n, alpha_ij, beta_im, z_i)

##### ρ ∈ M x N ####
cal_rho_MN <- function(t_m, T_star, A_m, mu_mi, Q_i) {
  rho_mn <- matrix(nrow = M,ncol = N)
  
  for (n in 1:N) {
    for (m in 1:M) {
      if (t_m[m] >= T_star[n]) {
        int_L_t <- as.interval(matrix(rep(c(A_m[m,1], T_star[n]), I), ncol = 2, byrow = T), L_open = F, R_open = T)
        rho_mn[m,n] <- sum(
          mu_mi[m,1:I]*is_subset(Q_i, int_L_t)
        )
      } else
      {
        rho_mn[m,n] <- 0
      }
    }
  }  
  rho_mn
}

# rho_mn <- cal_rho_MN(t_m, T_star, A_m, mu_mi, Q_i)


##### π ∈ U x N ####

# NEED FIX
cal_pi_UN <- function(t_u, T_star, A_u, eta_ui, Q_i) {
  pi_un <- matrix(nrow = U,ncol = N)
  
  for (n in 1:N) {
    for (u in 1:U) {
      if (t_u[u] >= T_star[n]) {
        int_L_t <- as.interval(matrix(rep(c(A_u[u,1], T_star[n]), I), ncol = 2, byrow = T), L_open = F, R_open = T)
        pi_un[u,n] <- sum(
           eta_ui[u,1:I]*is_subset(Q_i, int_L_t)
        )
      } else
      {
        pi_un[u,n] <- 0
      }
      
    }
  }  
  pi_un
}

# pi_un <- cal_pi_UN(t_u, T_star, A_u, eta_ui, Q_i)
##### σ ∈ C x N ####
cal_sigma_CN <- function(t_c, T_star, A_c, gamma_ci, Q_i) {
  sigma_cn <- matrix(nrow = C, ncol = N)
  
  for (n in 1:N) {
    for (c in 1:C) {
      if (t_c[c] >= T_star[n]) {
        int_L_t <- as.interval(matrix(rep(c(A_c[c,1], T_star[n]), I), ncol = 2, byrow = T), L_open = F, R_open = T)
        sigma_cn[c,n] <- sum(
          gamma_ci[c,1:I]*is_subset(Q_i, int_L_t)
        )
      } else
      {
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
e_24 <- function(c_k, mu_bar_ji, eta_ui, gamma_ci, N_star) {
  stopifnot(is.numeric(c_k), is.matrix(mu_bar_ji),
            is.matrix(eta_ui), is.matrix(gamma_ci),
            length(N_star) == 1, is.finite(N_star))
  
  I_diff <- I_mark - I
  mu_bar_ji <- mu_bar_ji[,(I+1):I_mark]
  eta_ui <- eta_ui[,(I+1):I_mark]
  gamma_ci <- gamma_ci[,(I+1):I_mark]
  
  if (!all(length(c_k) == I_diff, ncol(mu_bar_ji) == I_diff, ncol(eta_ui) == I_diff, ncol(gamma_ci) == I_diff)) {
    stop("mu_bar_ji, eta_ui, gamma_ci must each have K columns matching length(c_vec).")
  }
  
  (c_k + colSums(mu_bar_ji) + colSums(eta_ui) + colSums(gamma_ci)) / N_star
}

# λ_n from (25):
# λ_n = [ I(n ≤ N1) d_n + Σ_u I(t_{M+u} = t*_n) Σ_i η_{ui} ] /
#       [ Σ_m ρ_{mn} + Σ_u π_{un} + Σ_c σ_{cn} ]

e_25 <- function(d_n, t_u, T_star, eta_ui, rho_mn, pi_un, sigma_cn) {
  stopifnot(is.numeric(d_n), length(d_n) > 0)
  #N <- length(d_n)
  stopifnot(length(T_star) == N,
            ncol(rho_mn) == N,
            ncol(pi_un) == N,
            ncol(sigma_cn) == N,
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



