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
step_cdf <- function(times, masses, s_eval) {
  oo <- order(times, masses)
  times  <- times[oo]; masses <- masses[oo]
  cs <- cumsum(masses)
  # for each s, include all jumps with time ≤ s
  idx <- findInterval(s_eval, times, left.open = FALSE, rightmost.closed = TRUE)
  out <- ifelse(idx == 0, 0, cs[idx])
  pmin(out, 1)  # cap at 1, defensively
}

# helper: F(s-) at arbitrary t (strictly before t)
F_minus_at <- function(times_all, masses_all, t) {
  oo <- order(times_all, masses_all)
  ta <- times_all[oo]; ma <- masses_all[oo]
  cma <- cumsum(ma)
  idx <- findInterval(t, ta, left.open = TRUE, rightmost.closed = TRUE) # strict <
  ifelse(idx == 0, 0, pmin(cma[idx], 1))
}

# Main calculator
calc_F_and_hazards <- function(z_head, r, z_tail = numeric(0), Q = numeric(0),
                               s_eval, lambda = numeric(0), tstar = numeric(0),
                               l = NULL) {
  z_head <- as.numeric(z_head); r <- as.numeric(r)
  z_tail <- as.numeric(z_tail); Q <- as.numeric(Q)
  s_eval <- as.numeric(s_eval); lambda <- as.numeric(lambda); tstar <- as.numeric(tstar)
  
  stopifnot(length(z_head) == length(r))
  stopifnot(length(z_tail) == length(Q))
  
  # F12, F13, F on s_eval
  F12 <- if (length(r)) step_cdf(r, z_head, s_eval) else rep(0, length(s_eval))
  F13 <- if (length(Q)) step_cdf(Q, z_tail, s_eval) else rep(0, length(s_eval))
  Ftot <- pmin(F12 + F13, 1)
  
  # Build full jump lists for computing F(s-)
  all_times  <- c(r, Q)
  all_masses <- c(z_head, z_tail)
  
  # A12(s): denominators need F(l_i-) for each i
  if (is.null(l)) {
    # default: take l_i as an instant just before r_i; then F(l_i-) counts all jumps with time < r_i
    l_use <- r
  } else {
    l_use <- as.numeric(l)
    stopifnot(length(l_use) == length(r))
  }
  
  denom12 <- 1 - F_minus_at(all_times, all_masses, l_use)
  term12  <- ifelse(denom12 > 0, z_head / denom12, 0)
  
  # cumulative hazard A12 at s_eval: sum terms with r_i ≤ s
  A12 <- if (length(r)) {
    # cumulative sums aligned to ordered r
    oo <- order(r)
    r_o <- r[oo]; t_o <- term12[oo]
    cs  <- cumsum(t_o)
    idx <- findInterval(s_eval, r_o, left.open = FALSE, rightmost.closed = TRUE)
    ifelse(idx == 0, 0, cs[idx])
  } else rep(0, length(s_eval))
  
  # Λ13(s): denominators 1 - F(Q_k -)
  if (length(Q)) {
    denom13 <- 1 - F_minus_at(all_times, all_masses, Q)
    term13  <- ifelse(denom13 > 0, z_tail / denom13, 0)
    ooq <- order(Q)
    Q_o <- Q[ooq]; t13_o <- term13[ooq]
    cs13 <- cumsum(t13_o)
    idxq <- findInterval(s_eval, Q_o, left.open = FALSE, rightmost.closed = TRUE)
    A13 <- ifelse(idxq == 0, 0, cs13[idxq])
  } else {
    A13 <- rep(0, length(s_eval))
  }
  
  # Λ23(t) on t_eval = s_eval if you want; else evaluate at unique sorted union of tstar and s_eval
  t_eval <- sort(unique(s_eval))
  A23 <- if (length(tstar)) step_cdf(tstar, lambda, t_eval) else rep(0, length(t_eval))
  
  list(
    s_eval = s_eval,
    F12 = F12,
    F13 = F13,
    F = Ftot,
    Lambda12 = A12,
    Lambda13 = A13,
    t_eval_for_23 = t_eval,
    Lambda23 = A23
  )
}

