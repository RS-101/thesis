debugSource("aux_frydman.R")

#### Data generation ####
##### M : 1 -> 2 @ L,R -> 3 or N_tilde : ? @ t_m  ####
M <- 7 # interval transitions 1 -> 2.
N_tilde <- 3 # subset of 1:M 1 -> 2 -> 3 the rest of M is censored at 2.
if(N_tilde > M) stop("N_tilde should be lower than M")

id_m = seq(M)
in_N_tilde = id_m %in% sample(id_m, size = N_tilde)

L_m <- floor(runif(M, min = 1, max = 10))
R_m <- L_m + floor(runif(M, min = 1, max = 5))

t_m <- R_m+floor(runif(M, min = 0, max = 5))


##### K_tilde : 1 -> 3 @ e_k ####
K_tilde <- 6 # 1 -> 3
e_k <- floor(runif(K_tilde, min = 1, max = 10))

##### J : 1 -> ? @ s_j ####
J <- 5
s_j <- floor(runif(J, min = 1, max = 10))

##### U: 1 @ L_u- -> ? -> 3 @ t_u  ####
U <- 5 
L_u <- floor(runif(U, min = 1, max = 5)) # last time in state 1
t_u <- L_u + floor(runif(U, min = 2, max = 10))

##### C: 1 @ L_c- -> ? -> ? @ t_c  ####
C <- 5 
L_c <- floor(runif(C, min = 1, max = 5)) # last time in state 1
t_c <- L_u + floor(runif(C, min = 2, max = 10)) # censored

##### K: E* - Obs and potential 1 -> 3 ####
E_star <- unique(c(e_k,t_u))
c_k <- as.numeric(table(factor(c(e_k,t_u), levels = E_star)))
K <- length(E_star)

##### N: T* - Obs and potential 1 -> 2 -> 3 ####
T_star <- unique(c(t_m[in_N_tilde], t_u))
d_n <- as.numeric(table(factor(c(t_m[in_N_tilde], t_u), levels = T_star)))
N <- length(T_star)

##### Total: N* = M + U + C + K_tilde + J #### 
N_star <- M + U + C + K_tilde + J # Total count

##### Max 1 -> 2: M' = M + U + C ####
M_mark <- M + U + C # Max number through 2.

#### Creation of A sets ####

##### M: A_m := [L_m, R_m] ####
A_m <- as.interval(matrix(c(L_m, R_m), ncol = 2, byrow = F))

##### W: M < m <= W := M + U, R_{M+u} = t_{M+u} ####
A_u <- as.interval(matrix(c(L_u,t_u), ncol = 2, byrow = F))

##### M': W := M + U < m <= M', R_{W+c} = t_{W+c} ####
A_c <- as.interval(matrix(c(L_c,t_c), ncol = 2, byrow = F))

##### full_A_m: A_m ∪ A_u ∪ A_c ####
full_A_m <- as.interval(rbind(A_m,A_u,A_c))

##### A := ⋃_{m=1}^{M'} A_m ####
A_union <- get_interval(full_A_m)
A_union <- get_interval(as.interval(matrix(c(1,4),nrow = 1,byrow = T)))

#### Data manipulation ####
##### I: Q_i = [l_i,r_i] ####

# s_max = max(s_j, 1 <= j <= J)
s_max = max(s_j)

# R_max = max(R_m, 1 <= m <= W)
R_max = max(A_m[,2], A_u[,2])

# e*_max = max(e*_k, 1 <= k <= K)
e_star_max = max(E_star)

# L_bar ={L_m, 1 <= m <= M'} ∪ {T* ∩ A} ∪ {S_J ∩ A} ∪ {s_max : s_max > R_max ∨ e*_max}
L_bar <- c(full_A_m[,1],
           intersect(A_union,T_star),
           intersect(s_j, A_union),
           na.omit(ifelse(s_max > max(R_max,e_star_max), s_max, NA)))

# R_bar = {R_m, 1 <= m <= W} ∪ {∞}
R_bar = c(full_A_m[1:(M+U),2], Inf)

L_bar_sort_distinct <- sort(unique(L_bar))
R_bar_sort_distinct <- sort(unique(R_bar))

# !!!DANGER I AM UNSURE ABOUT THE CREATION OF Q!!!
Q_i <- make_Q(L_bar_sort_distinct, R_bar_sort_distinct)

Q <- get_interval(Q_i)

I <- nrow(Q)


##### I'= K + I: Q_i' = e*_i-I ####
Q_i_mark <- E_star

##### Q_full = list ####
Q_full <- list(Q, Q_i_mark)

##### C: s_J+c = t_W+c ####
s_j_c <- t_c

s_j_full <- c(s_j,s_j_c)


##### N: lambda_n and I': z_i ####
# Comment: I believe we have I' z_i's and N: lambda_n
I_mark <- I + K

z_i_0 <- numeric(I_mark)
lambda_n_0 <- numeric(N)

#### Creating functions from (18) - (25) ####
# !!!DANGER SHOULD IT BE SHARP OR NOT!!!
# alpha is a I' x (J+C) matrix ~ alpha_ij = I(Q_i_full subset [s_j_full,inf)).

alpha <- function(Q_full, s_j_full) {
  Q_i <- Q_full[[1]]
  Q_i_mark <- Q_full[[2]]
  alpha <- matrix(numeric(I_mark*(J+C)), nrow = I_mark, ncol = J+C)
  for(i in 1:I_mark) {
    if(i <= I) {
      alpha[i,] <- s_j_full <= Q_i[i,1]
    } else {
      alpha[i,] <- s_j_full <= Q_i_mark[i-I]
    }
  }
  alpha
}

alpha_ij <- alpha(Q_full, s_j_full)

# !!!DANGER SHOULD IT BE SHARP OR NOT!!!
# beta is a I x M' matrix ~ beta_im = I(Q_i subset A_m).
beta <- function(Q_i, full_A_m) {
  beta <- matrix(numeric(I*M_mark), nrow = I, ncol = M_mark)
  for(i in 1:I) {
    beta[i,] <- full_A_m[,1] <= Q_i[i,1] & Q_i[i,2] <= full_A_m[,2]
  }
  beta
}
beta_im <- beta(Q_i, full_A_m)




##### μ_mi(z,λ,β,Q,A_m) ∈ M x I, ####  
# mu_m is a matrix M x I
mu <- function(z, lambda, beta_im, Q_i, A_m, T_star) {
  beta_im <- beta_im[1:I, 1:M]
  if(!all(dim(beta_im) == c(I,M))) stop("beta wrong dim")
  
  mu <- matrix(numeric(I*M), nrow = M, ncol = I)
  
  r <- Q[,2]
  R <- A_m[,2]
  
  for (m in 1:M) {
   prod_res <- unlist(lapply(r, function(r_i) product_over_t_stars(intervals = as.interval(
      matrix(c(r_i,R[m]),ncol = 2), 
      L_open = T, 
      R_open = F),
      T_star = T_star,
      lambda = lambda)))
    denum <- sum(beta_im[,m]*z[1:I]*prod_res)
    for (i in 1:I) {
      num <- beta_im[i,m]*z[i]*prod_res[i]
      mu[m,i] <- ifelse(num != 0, num/denum, 0)
    }
  }
  mu
}
mu(z, lambda, beta_im, Q_i, A_m, T_star )

mu_bar <- function() {
  
}

neta <- function() {
  
}

gamma <- function() {
  
}


# equation (23) i <= I
e_23 <- function(z,lambda) {
  
  # we get a vector of I values
  
  # bar_mu_j is a matrix J x I
  # neta_u is a matrix U x I 
  # gamma_c is a matrix C x I
  num <- sum(mu_m) + sum(mu_bar_j) + sum(neta_u) + sum(gamma_c)
  denum <- N_star
  
  return(z)
} 

# equation (24) I < i <= I'
e_24 <- function(z, lambda) {
  return(z)
}

e_25 <- function(z,lambda) {
  
  return(lambda)
}
