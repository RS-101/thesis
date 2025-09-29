# Compute F̂12(s), F̂13(s), F̂(s)=F̂12+F̂13 and cumulative hazards Λ̂12(s), Λ̂13(s),
# plus Λ̂23(t) from (text above) given (ẑ, λ̂) and jump times.
#
# Inputs
#   z_head : length I,   ẑ_1..ẑ_I        (jumps tied to r_i)
#   r      : length I,   r_1..r_I          (jump times for F12)
#   z_tail : length K,   ẑ_{I+1}..ẑ_{I+K} (jumps tied to Q)
#   Q      : length K,   Q_{I+1}..Q_{I+K}  (jump times for F13)
#   s_eval : vector of s where to evaluate F and Λ for states 1→2 and 1→3
#   lambda : length N, λ̂_1..λ̂_N
#   tstar  : length N, t*_1..t*_N          (jump times for Λ23)
#
# Notes
# - F̂12(s) = Σ_{i: r_i ≤ s} ẑ_i
# - F̂13(s) = Σ_{k: Q_k ≤ s} ẑ_{I+k}
# - F̂(s)   = F̂12(s) + F̂13(s)
# - Λ̂12(s) = Σ_{i: r_i ≤ s} ẑ_i / {1 - F̂(l_i-)}   with l_i- taken as F̂ just before l_i.
#            If l_i are not available, we take l_i = r_{i-1}^+ so F̂(l_i-) = Σ_{j<i} ẑ_j + Σ_{Q<Q_i} ẑ_Q.
#            You can pass explicit l if you have them.
# - Λ̂13(s) = Σ_{k: Q_k ≤ s} ẑ_{I+k} / {1 - F̂(Q_k-)}
# - Λ̂23(t) = Σ_{n: t*_n ≤ t} λ̂_n
#
# If you have explicit l (the left-interval points), pass them; else it will
# approximate l_i- by the time just before r_i using the rule above.

na0 <- function(x) { x[is.na(x)] <- 0; x }

# helper: right-continuous step CDF from (times, masses) at s_eval
step_cdf <- function(grid_points, times, masses, ...) {
  time_order <- order(times)
  times <- times[time_order]
  masses <- masses[time_order]
  cs <- cumsum(masses)
  idx <- findInterval(grid_points, times, ...)
  c(0,cs)[idx+1]
}

# Main calculator
calc_F_and_hazards <- function(grid_points, z_i, lambda_n, Q_i, Q_i_mark, T_star, E_star) {
  
  
  I <- nrow(Q_i)
  I_mark <- I + length(Q_i_mark)
  # F12, F13, F on s_eval
  
  # use intercept for F12 instead of just right or left endpoint
  
  F12 <- step_cdf(grid_points, times = Q_i[1:I, 2], masses = z_i[1:I])
  F13 <- step_cdf(grid_points, times = Q_i_mark, masses = z_i[(I+1):I_mark])
  F_total <- F12 + F13
  
  A23 <- step_cdf(grid_points, T_star, masses = lambda_n)
  
  
  F12_at_l_i <- step_cdf(Q_i[,1]-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I])
  # A12(s): denominators need F(l_i-) for each i
  denom12 <- 1 - F12_at_l_i
  term12  <- ifelse(denom12 > 0, z_i[1:I] / denom12, 0)
  A12 <- step_cdf(grid_points, times = Q_i[,2], term12)
  
  F_total_at_e_k <- step_cdf(E_star-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I]) + 
    step_cdf(E_star-1e-6, times = Q_i_mark, masses = z_i[(I+1):I_mark])
  # A12(s): denominators need F(l_i-) for each i
  denom13 <- 1 - F_total_at_e_k
  term13  <- ifelse(denom12 > 0, z_i[(I+1):I_mark] / denom13, 0)
  A13 <- step_cdf(grid_points, times = E_star, term13)
  
  list(
    grid_points = grid_points,
    F12 = F12,
    F13 = F13,
    F = F_total,
    Lambda12 = A12,
    Lambda13 = A13,
    Lambda23 = A23
  )
}

