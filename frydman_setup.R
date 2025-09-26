
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

z_i <- runif(I_mark)
lambda_n <- runif(N)
