# ggplot2 rewrite (uses ggplot2 + patchwork)
library(ggplot2)
library(dplyr)
library(patchwork)
library(data.table)

#### Joly data to Frydman data
setup_data_to_list_format <- function(data, add_r_format = F) {

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

  get_interval <- function(m, L_open = F, R_open = F) {
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





  
  to_mat <- function(x) if (is.matrix(x)) x else as.matrix(unclass(x))
  
  data <- as.data.table(data)
  if(isFALSE("id" %in% names(data))) {data[, id:= .I]}
  
  stopifnot(all(sort(names(data)) == sort(c("id", "V_0", "V_healthy",
                                            "V_ill", "T_obs", "status"))))
  # case 1
  case_1 <- data[as.integer(status) == 1]
  case_1_exact <- case_1[V_healthy == T_obs]
  
  J <- case_1_exact[,.N]
  s_j <- case_1_exact$T_obs
  
  
  case_1_rest <- case_1[!(V_healthy == T_obs)]
  C <- case_1_rest[,.N]
  L_c <- case_1_rest$V_healthy
  t_c <- case_1_rest$T_obs
  
  # case 2
  case_2 <- data[as.integer(status) == 2]
  case_2_exact <- case_2[V_healthy == T_obs]
  K_tilde <- case_2_exact[,.N]
  e_k <- case_2_exact$T_obs
  
  case_2_rest <- case_2[!(V_healthy == T_obs)]
  U <- case_2_rest[,.N]  
  L_u <- case_2_rest$V_healthy
  t_u <- case_2_rest$T_obs
  
  # case 4
  case_4 <- data[as.integer(status) == 4]
  N_tilde <- case_4[,.N]  
  t_m_in_N_tilde <- case_4$T_obs
  
  # case 3
  case_3 <- data[as.integer(status) == 3]
  case_3_4 <- data[as.integer(status) %in% c(3,4)]
  M <- case_3_4[,.N]
  
  L_m <- case_3_4$V_healthy
  R_m <- case_3_4$V_ill
  t_m <- case_3_4$T_obs
  
  ##### K: E* - Obs and potential 1 -> 3 ####
  E_star <- unique(c(e_k, t_u))
  #c_k <- as.numeric(table(factor(c(e_k, t_u), levels = E_star)))
  # DANGER CHAT SAYS THAT THIS IS CORRECT
  # as sum(c_k) then is equal to K_tilde, otherwise K < sum(c_k) <= K + U, with = when e_k ∩ t_u = ∅
  c_k <- as.numeric(table(factor(c(e_k, t_u), levels = E_star)))
  K <- length(E_star)
  
  ##### N: T* - Obs and potential entry to state 3 from state 2: 1 -> 2 -> 3 ####
  t_star_n <- unique(c(t_m_in_N_tilde, t_u))
  d_n <- as.numeric(table(factor(c(t_m_in_N_tilde, t_u), levels = t_star_n)))
  
  N1_obs_of_T_star <- length(unique(t_m_in_N_tilde))
  U_pos_obs_of_T_star <- length(setdiff(unique(t_u), unique(t_m_in_N_tilde)))
  
  N <- length(t_star_n)
  
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
  s_max <- ifelse(length(s_j) > 0 ,max(s_j),0)
  
  # R_max = max(R_m, 1 <= m <= W)
  R_max <- max(A_m[, 2], A_u[, 2])
  
  # e*_max = max(e*_k, 1 <= k <= K)
  e_star_max <- max(E_star)
  
  # L_bar ={L_m, 1 <= m <= M'} ∪ {T* ∩ A} ∪ {S_J ∩ A} ∪ {s_max : s_max > R_max ∨ e*_max}
  L_bar <- c(
    full_A_m[, 1],
    intersect.interval(A_union, t_star_n),
    intersect.interval(A_union, s_j),
    na.omit(ifelse(s_max > max(R_max, e_star_max), s_max, NA))
  )
  
  # R_bar = {R_m, 1 <= m <= W} ∪ {∞}
  R_bar <- c(full_A_m[1:(M + U), 2], Inf)
  
  # !!!DANGER I AM UNSURE ABOUT THE CREATION OF Q!!!
  Q_i <- make_Q(L_bar, R_bar)
  
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
  
  data_list <- list(
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
    E_star = E_star, t_star_n = t_star_n, c_k = c_k, d_n = d_n,
    L_bar = L_bar, R_bar = R_bar, s_j_c = s_j_c, s_j_full = s_j_full,
    Q = Q, Q_i_mark = Q_i_mark, A_union = A_union,
    
    # 2-col “interval” matrices
    A_m = to_mat(A_m), A_u = to_mat(A_u), A_c = to_mat(A_c),
    full_A_m = to_mat(full_A_m), Q_i = to_mat(Q_i)
  )
  
  if(add_r_format) {
    return(list(cpp_data = data_list, 
         r_data = list(
           Q_full   = Q_full,
           s_j_full = s_j_full, 
           Q_i      = Q_i,
           full_A_m = full_A_m, 
           A_m      = A_m,
           A_u      = A_u,
           A_c      = A_c, 
           t_star_n   = t_star_n,
           E_star   = E_star,
           t_m      = t_m,
           t_u      = t_u,
           t_c      = t_c,
           N_star   = N_star,
           d_n      = d_n,
           c_k      = c_k,
           I_mark   = I_mark,
           J        = J, 
           M        = M,
           W        = W,
           K_tilde  = K_tilde,
           N1_obs_of_T_star = N1_obs_of_T_star
         )))
  }
  
  data_list
}
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
  cs <- pmax(cumsum(masses),0)
  idx <- findInterval(grid_points, times, ...)
  c(0,cs)[idx+1]
}

# Main calculator
calc_F_and_hazards <- function(grid_points, z_i, lambda_n, Q_i, Q_i_mark, t_star_n, E_star, as_function = FALSE) {
  
  
  I <- nrow(Q_i)
  I_mark <- I + length(Q_i_mark)
  # F12, F13, F on s_eval
  
  # use intercept for F12 instead of just right or left endpoint
  
  F12 <- step_cdf(grid_points, times = Q_i[1:I, 2], masses = z_i[1:I])
  F13 <- step_cdf(grid_points, times = Q_i_mark, masses = z_i[(I+1):I_mark])
  F_total <- F12 + F13
  
  A23 <- step_cdf(grid_points, t_star_n, masses = lambda_n)
  
  
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

  
  # Turn vectors into step functions over grid_points
  stepify <- function(x, y, side = c("left", "right"), extend = TRUE) {
    side <- match.arg(side)
    f <- if (side == "left") 0 else 1
    rule <- if (extend) 2 else 1  # 2 = hold ends constant, 1 = NA outside
    approxfun(x, y, method = "constant", f = f, rule = rule)
  }
  
  list(
    as_points = list(
      F12       = stepify(grid_points, F12, side = "left"),
      F13       = stepify(grid_points, F13, side = "left"),
      F         = stepify(grid_points, F_total, side = "left"),
      Lambda12  = stepify(grid_points, A12, side = "left"),
      Lambda13  = stepify(grid_points, A13, side = "left"),
      Lambda23  = stepify(grid_points, A23, side = "left")
    ),
    as_functions = list(
        grid_points = grid_points,
        F12 = F12,
        F13 = F13,
        F = F_total,
        Lambda12 = A12,
        Lambda13 = A13,
        Lambda23 = A23
      ))
}


plot_estimators_gg <- function(res) {
  # Common theme
  thm <- theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_line(),
          plot.title = element_text(face = "bold"))
  
  # F12
  df_F12 <- data.frame(s = res$grid_points, y = res$F12) |> arrange(s)
  p1 <- ggplot(df_F12, aes(x = s, y = y)) +
    geom_step(linewidth = 1, color = "blue") +
    labs(x = "s", y = expression(hat(F)[12](s)), title = "Estimator F12(s)") +
    thm
  
  # F13
  df_F13 <- data.frame(s = res$grid_points, y = res$F13) |> arrange(s)
  p2 <- ggplot(df_F13, aes(x = s, y = y)) +
    geom_step(linewidth = 1, color = "red") +
    labs(x = "s", y = expression(hat(F)[13](s)), title = "Estimator F13(s)") +
    thm
  
  # F
  df_F <- data.frame(s = res$grid_points, y = res$F) |> arrange(s)
  p3 <- ggplot(df_F, aes(x = s, y = y)) +
    geom_step(linewidth = 1, color = "purple") +
    labs(x = "s", y = expression(hat(F)(s)), title = "Estimator F(s)") +
    thm
  
  # Hazards Λ12, Λ13, Λ23
  df_haz <- dplyr::bind_rows(
    data.frame(s = res$grid_points,         y = res$Lambda12, which = "L12"),
    data.frame(s = res$grid_points,         y = res$Lambda13, which = "L13"),
    data.frame(s = res$grid_points,         y = res$Lambda23, which = "L23")
  ) |> arrange(s)
  
  p4 <- ggplot(df_haz, aes(x = s, y = y, color = which)) +
    geom_step(linewidth = 1, alpha = 0.5) +
    scale_color_manual(
      values = c(L12 = "blue", L13 = "red", L23 = "darkgreen"),
      labels = c(
        L12 = expression(hat(Lambda)[12]),
        L13 = expression(hat(Lambda)[13]),
        L23 = expression(hat(Lambda)[23])
      )
    ) +
    labs(x = "s", y = "Cumulative hazards", color = NULL, title = "Λ estimators") +
    thm
  
  # Arrange 2x2 grid
  (p1 | p2) / (p3 | p4)
}

#' Plot estimated cumulative hazards (Lambda12/13/23) against true cumulative hazards
#'
#' @param estimators A list with elements:
#'   - grid_points: numeric vector of evaluation points (increasing)
#'   - Lambda12, Lambda13, Lambda23: numeric vectors (same length as grid_points)
#' @param true_hazards A list (e.g. sim_dat$true_data_generation$hazards) with
#'   functions a12(t), a13(t), a23(s) returning hazards.
#' @param xlim Optional x-axis limits. Defaults to range(estimators$grid_points).
#' @param lower Optional lower limit for cumulative integration of true hazards.
#'   Defaults to min(xlim).
#' @param title Plot title.
#' @param ... Passed to ggplot2::geom_line (e.g. linewidth).
#'
#' @return A ggplot object.
plot_cumhazards_est_vs_true <- function(estimators,
                                        true_hazards,
                                        xlim = range(estimators$grid_points, na.rm = TRUE),
                                        lower = min(xlim),
                                        title = "Cumulative hazard: estimate vs true",
                                        ...) {
  stopifnot(is.list(estimators),
            is.numeric(estimators$grid_points),
            all(c("Lambda12", "Lambda13", "Lambda23") %in% names(estimators)),
            is.list(true_hazards))
  
  # Extract grid and estimated cumulative hazards --------------------------------
  grid <- as.numeric(estimators$grid_points)
  ord  <- order(grid)
  grid <- grid[ord]
  
  est_df <- data.frame(
    x   = rep(grid, 3L),
    y   = c(estimators$Lambda12[ord],
            estimators$Lambda13[ord],
            estimators$Lambda23[ord]),
    fun = rep(c("a12", "a13", "a23"), each = length(grid)),
    type = "estimate",
    stringsAsFactors = FALSE
  )
  
  # Helpers ----------------------------------------------------------------------
  integrate_on_grid <- function(f, grid, lower) {
    vapply(
      grid,
      function(t1) {
        if (t1 <= lower) return(0)
        stats::integrate(function(u) f(u), lower = lower, upper = t1, rel.tol = 1e-6)$value
      },
      numeric(1)
    )
  }
  
  # Build true cumulative hazards by integrating the true hazards ----------------
  # Support both true_hazards$hazards$list_of_functions and the list itself
  as_fun_list <- function(obj) {
    f <- if (!is.null(obj$hazards)) obj$hazards else obj
    f[unlist(lapply(f, is.function))]
  }
  th <- as_fun_list(true_hazards)
  
  req <- c("a12", "a13", "a23")
  if (!all(req %in% names(th))) {
    stop("true_hazards must contain functions named a12, a13, a23.")
  }
  
  true_df <- do.call(
    rbind,
    lapply(req, function(nm) {
      data.frame(
        x   = grid,
        y   = integrate_on_grid(th[[nm]], grid, lower = lower),
        fun = nm,
        type = "true",
        stringsAsFactors = FALSE
      )
    })
  )
  
  df <- rbind(est_df, true_df)
  df$fun <- factor(df$fun, levels = req, labels = req)
  
  # Plot -------------------------------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = x, y = y, color = fun, linetype = type)
  ) +
    ggplot2::geom_line(alpha = 0.9, ...) +
    ggplot2::scale_linetype_manual(values = c(estimate = "dashed", true = "solid")) +
    ggplot2::labs(
      title = title,
      x = "time",
      y = "cumulative hazard",
      color = "transition",
      linetype = NULL
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::theme_minimal()
  
  p
}

# Example (assuming objects exist):
# p <- plot_cumhazards_est_vs_true(estimators, sim_dat$true_data_generation$hazards)
# print(p)



