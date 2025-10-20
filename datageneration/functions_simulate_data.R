# --- Verification helper for the illness–death model -------------------------
# Requires: simulate_illness_death(a12, a13, a23, ...) from earlier.
# What it does:
#   1) Simulate n trajectories.
#   2) Compute theoretical PDFs/CDFs:
#        - First event time:        f_first(t) = S(t) * (a12(t)+a13(t)),  F_first(t) = 1 - S(t)
#        - Illness time (given illness occurs):  f12|ill(t) ∝ S(t)*a12(t)
#        - Direct death time (given direct death): f13|dir(t) ∝ S(t)*a13(t)
#        - Ill→death elapsed time s:
#             * if a23 = a23(s): f23(s) = S23(s) * a23(s), F23(s) = 1 - S23(s)
#             * if a23 = a23(t_abs, t_entry): mixture over observed entry times
#   3) Plot histograms + overlays and return functions for the CDF/PDFs.


# Simulate an illness–death model with arbitrary hazards
# - Everyone starts in state 1 (healthy)
# - Independently draw T12 ~ a12(t) and T13 ~ a13(t); the minimum determines the first transition
# - If T12 <= T13, enter state 2 at time T12 and then draw time to death from illness
#   conditionally on entry, using hazard a23.
#
# a23 can be specified in two ways:
#   (A) a function of elapsed time since illness:        a23 <- function(s) {...}
#   (B) a function of absolute time and entry time:      a23 <- function(t_abs, t_entry) {...}
# The simulator detects which form you provided.
#
# Returns a data.frame with id, entry2 (time of illness or Inf),
# t13_direct (time of direct death from state 1), t23_after (time from illness to death),
# death_time (absolute), and path ("1->3" or "1->2->3").


library(ggplot2)
library(patchwork)

simulate_illness_death <- function(
    n,
    a12, a13, a23,      # all are functions of ABSOLUTE calendar time t
    t0 = 0,
    tmax = Inf,          # if finite: events after tmax are treated as never occurring (Inf)
    init_step = 1,       # initial bracket step for inversion search
    int_abs_tol = 1e-8,  # absolute tolerance for integrate()
    root_tol = 1e-8,     # tolerance for uniroot()
    max_doublings = 60   # safety cap for growing the search bracket
) {
  
  # --- Safeguards -------------------------------------------------------------
  if (length(formals(a12)) != 1L || length(formals(a13)) != 1L)
    stop("a12 and a13 must be functions of one argument: absolute time t.")
  if (length(formals(a23)) != 1L)
    stop("a23 must be a function of one argument: absolute time t (calendar-time hazard).")
  
  # --- Helper: cumulative hazard from t_start to t_end ------------------------
  cumhaz <- function(h, t_start, t_end) {
    if (t_end <= t_start) return(0)
    res <- stats::integrate(function(u) h(u), lower = t_start, upper = t_end,
                            stop.on.error = TRUE, abs.tol = int_abs_tol)
    as.numeric(res$value)
  }
  
  # --- Helper: draw absolute event time by inverting the cumulative hazard ----
  # Given hazard h(t) and starting clock t_start, draw T ≥ t_start with
  #  P(T > t | T ≥ t_start) = exp(-∫_{t_start}^t h(u) du).
  draw_time_from_hazard <- function(h, t_start) {
    target <- stats::rexp(1)  # -log(U) ~ Exp(1)
    
    # grow an upper bound ub until ∫ h >= target, or we hit tmax
    lb <- t_start
    step <- init_step
    ub <- lb + step
    ch <- cumhaz(h, lb, ub)
    
    n_doubles <- 0
    while (is.finite(ch) && ch < target && ub < tmax && n_doubles < max_doublings) {
      lb <- ub
      step <- step * 2
      ub <- lb + step
      ch <- ch + cumhaz(h, lb, ub)  # accumulate ∫ on the new segment
      n_doubles <- n_doubles + 1
    }
    
    # if we couldn’t reach target before tmax, the event never occurs
    if (!is.finite(ch) || (ub >= tmax && ch < target)) return(Inf)
    
    # root-find F(t) = ∫_{t_start}^{t} h(u) du - target = 0 on [lb0, ub0]
    F <- function(t) cumhaz(h, t_start, t) - target
    lb0 <- t_start
    ub0 <- ub
    # F(lb0) < 0 always (target > 0); ensure F(ub0) >= 0
    val_ub <- F(ub0)
    while (val_ub < 0 && ub0 < tmax && n_doubles < max_doublings + 20) {
      lb0 <- ub0
      ub0 <- ub0 + step
      val_ub <- F(ub0)
      n_doubles <- n_doubles + 1
    }
    if (val_ub < 0) return(Inf)
    
    uniroot(F, lower = lb0, upper = ub0, tol = root_tol)$root  # absolute time
  }
  
  # --- Main simulation loop ---------------------------------------------------
  out <- vector("list", n)
  for (i in seq_len(n)) {
    
    # Draw absolute times from state 1 (competing risks, calendar time)
    # Time-change theorem: if E ~ Exp(1), then T = inf{ t ≥ t0 : ∫_{t0}^t a_{1k}(u) du ≥ E }
    t12 <- draw_time_from_hazard(a12, t0)
    t13 <- draw_time_from_hazard(a13, t0)
    
    if (t12 <= t13) {
      # Illness occurs first at entry2 = t12
      entry2 <- t12
      
      # From state 2, with calendar-time hazard a23(t):
      # Conditional on entry at u, T | U=u has survival exp(-∫_u^T a23(v) dv).
      # We sample the ABSOLUTE death time:
      death_abs <- draw_time_from_hazard(a23, entry2)
      
      # Store absolute death time and the sojourn S = T - U
      if (is.finite(death_abs)) {
        t23_after <- death_abs - entry2
        death_time <- death_abs
      } else {
        t23_after <- NA_real_   # censored after entry2
        death_time <- Inf
      }
      
      out[[i]] <- list(
        id = i,
        entry2 = entry2,
        t13_direct = t13,       # drawn but not realized
        t23_after = t23_after,  # sojourn in state 2 (NA if censored)
        death_time = death_time,
        path = "1->2->3"
      )
      
    } else {
      # Direct death from state 1 at t13
      out[[i]] <- list(
        id = i,
        entry2 = Inf,
        t13_direct = t13,
        t23_after = NA_real_,
        death_time = t13,
        path = "1->3"
      )
    }
  }
  
  df <- do.call(rbind, lapply(out, as.data.frame))
  rownames(df) <- NULL
  df$path <- factor(df$path, levels = c("1->3", "1->2->3"))
  df
}

# Plots BOTH:
#  - sojourn density f_S(s) for Ill->Death elapsed time (mixture over entry times)
#  - calendar-time density f_T(t | path=1->2->3) for death time among those who went through illness
verify_illness_death <- function(
    n,
    a12, a13, a23,      # all: function(t_abs); a23 is CALENDAR-TIME hazard
    t0 = 0,
    t_grid = NULL,
    s_grid = NULL,
    ngrid = 600,
    seed = NULL,
    sim_fun = simulate_illness_death,
    print_plot = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  
  trapz_cum <- function(x, y) {
    n <- length(x)
    out <- numeric(n)
    for (i in 2:n) out[i] <- out[i-1] + 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    out
  }
  trapz_cols <- function(x, Y) {
    # integrate each column of Y over x by trapezoid
    dx <- diff(x)
    Yu <- Y[-1, , drop = FALSE]
    Yd <- Y[-nrow(Y), , drop = FALSE]
    colSums((Yu + Yd) * rep(dx, times = ncol(Y)) / 2)
  }
  
  # 1) Simulate data (histograms)
  sim <- sim_fun(n = n, a12 = a12, a13 = a13, a23 = a23, t0 = t0)
  
  first_time <- ifelse(sim$path == "1->2->3", sim$entry2, sim$t13_direct)
  got_ill   <- is.finite(sim$entry2)
  dir_death <- sim$path == "1->3"
  s_23      <- sim$t23_after[!is.na(sim$t23_after)]
  death_123 <- sim$death_time[sim$path == "1->2->3"]
  
  # 2) Grids
  if (is.null(t_grid)) {
    t_end <- max(first_time[is.finite(first_time)], na.rm = TRUE)
    t_end <- max(t_end, stats::quantile(first_time, 0.995, na.rm = TRUE))
    t_end <- max(t_end, stats::quantile(sim$death_time, 0.995, na.rm = TRUE))
    t_end <- t_end * 1.25 + 1e-8
    t_grid <- seq(t0, t_end, length.out = ngrid)
  }
  if (is.null(s_grid)) {
    if (length(s_23)) {
      s_end <- max(stats::quantile(s_23, 0.995, na.rm = TRUE), max(s_23, na.rm = TRUE))
      s_end <- s_end * 1.25 + 1e-8
    } else {
      s_end <- 1
    }
    s_grid <- seq(0, s_end, length.out = ngrid)
  }
  
  # 3) First-event theory on t_grid
  h_sum  <- a12(t_grid) + a13(t_grid)
  H_sum  <- trapz_cum(t_grid, h_sum)
  S1     <- exp(-H_sum)
  f12    <- S1 * a12(t_grid)
  f13    <- S1 * a13(t_grid)
  f_first<- f12 + f13
  F12    <- trapz_cum(t_grid, f12)
  F13    <- trapz_cum(t_grid, f13)
  F_first<- 1 - S1
  P12    <- F12[length(F12)]
  P13    <- F13[length(F13)]
  
  f12_cond <- if (P12 > 0) f12 / P12 else rep(0, length(f12))
  f13_cond <- if (P13 > 0) f13 / P13 else rep(0, length(f13))
  F12_cond <- if (P12 > 0) F12 / P12 else F12
  F13_cond <- if (P13 > 0) F13 / P13 else F13
  
  # 4) Calendar-time a23: precompute A23(t)
  tA_max   <- max(t_grid) + max(s_grid)
  tA_grid  <- seq(min(t_grid), tA_max, length.out = max(ngrid, length(t_grid)))
  a23_vals <- a23(tA_grid)
  A23      <- trapz_cum(tA_grid, a23_vals)
  A23_fun  <- stats::approxfun(tA_grid, A23, rule = 2)
  a23_fun  <- stats::approxfun(tA_grid, a23_vals, rule = 2)
  
  te <- t_grid
  s  <- s_grid
  
  # 4a) Sojourn s = t - te (mixture over entry-time density f12|ill(te))
  te_plus_s   <- outer(te, s, "+")
  A_teps      <- A23_fun(te_plus_s)
  A_te_mat    <- matrix(A23_fun(te), nrow = length(te), ncol = length(s), byrow = FALSE)
  S23_te_s    <- exp(-(A_teps - A_te_mat))                 # S23(s | te)
  a23_teps    <- a23_fun(te_plus_s)                        # a23(te + s)
  dens_te_s   <- matrix(f12_cond, nrow = length(te), ncol = length(s), byrow = FALSE)
  integrand_S <- S23_te_s * dens_te_s
  integrand_f <- S23_te_s * a23_teps * dens_te_s
  S23_mix     <- if (P12 > 0) trapz_cols(te, integrand_S) else rep(1, length(s))
  f23_mix     <- if (P12 > 0) trapz_cols(te, integrand_f) else rep(0, length(s))
  F23_mix     <- 1 - S23_mix
  
  # 4b) Death-time density conditional on path 1->2->3 (calendar time)
  # f_T(t | 1->2->3) = ∫_{u<=t} f12|ill(u) * S23(t | u) * a23(t) du
  A_t_col     <- matrix(A23_fun(t_grid), nrow = length(te), ncol = length(t_grid), byrow = TRUE)
  A_te_col    <- matrix(A23_fun(te),     nrow = length(te), ncol = length(t_grid), byrow = FALSE)
  S23_tu      <- exp(-(A_t_col - A_te_col))                # S23(t | u)
  mask_ut     <- outer(te, t_grid, function(u, t) as.numeric(u <= t))
  a23_t_col   <- matrix(a23_fun(t_grid), nrow = length(te), ncol = length(t_grid), byrow = TRUE)
  dens_te_t   <- matrix(f12_cond, nrow = length(te), ncol = length(t_grid), byrow = FALSE)
  integrand_ft<- S23_tu * a23_t_col * dens_te_t * mask_ut
  f_death_123 <- if (P12 > 0) trapz_cols(te, integrand_ft) else rep(0, length(t_grid))
  F_death_123 <- trapz_cum(t_grid, f_death_123)
  
  # 5) Plots
  
  # Helper for FD binwidth (like breaks = "FD")
  fd_binwidth <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NULL)
    bw <- 2 * IQR(x) / (length(x)^(1/3))
    if (is.finite(bw) && bw > 0) bw else NULL
  }
  
  # Theoretical curves as data frames
  df_A  <- data.frame(t = t_grid, f = f_first)
  df_B  <- if (exists("f12_cond")) data.frame(t = t_grid, f = f12_cond) else NULL
  df_C  <- if (exists("f13_cond")) data.frame(t = t_grid, f = f13_cond) else NULL
  df_D  <- data.frame(s = s_grid,  f = f23_mix)
  df_E  <- data.frame(t = t_grid,  f = f_death_123)
  
  # (A) First event time
  pA <- ggplot() +
    geom_histogram(aes(x = first_time, y = after_stat(density)),
                   binwidth = fd_binwidth(first_time),
                   color = "grey40", fill = "grey80") +
    geom_line(data = df_A, aes(t, f, color = "theory f_first(t)"), linewidth = 1) +
    labs(title = "First event time (min of illness/death)", x = "t", y = "Density", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  # (B) Illness time | illness occurs
  if (any(got_ill)) {
    xB <- sim$entry2[got_ill]
    pB <- ggplot() +
      geom_histogram(aes(x = xB, y = after_stat(density)),
                     binwidth = fd_binwidth(xB),
                     color = "grey40", fill = "grey80") +
      geom_line(data = df_B, aes(t, f, color = "theory f12|ill(t)"), linewidth = 1) +
      labs(title = "Illness time | illness occurs", x = "t", y = "Density", color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  } else {
    pB <- ggplot() +
      xlim(0, 1) + ylim(0, 1) + theme_void(base_size = 11) +
      geom_text(aes(0.5, 0.5, label = "No illness events observed"))
  }
  
  # (C) Direct death time | direct death
  if (any(dir_death)) {
    xC <- sim$death_time[dir_death]
    pC <- ggplot() +
      geom_histogram(aes(x = xC, y = after_stat(density)),
                     binwidth = fd_binwidth(xC),
                     color = "grey40", fill = "grey80") +
      geom_line(data = df_C, aes(t, f, color = "theory f13|dir(t)"), linewidth = 1) +
      labs(title = "Direct death time | direct death", x = "t", y = "Density", color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  } else {
    pC <- ggplot() +
      xlim(0, 1) + ylim(0, 1) + theme_void(base_size = 11) +
      geom_text(aes(0.5, 0.5, label = "No direct deaths observed"))
  }
  
  # (D) Ill → death elapsed time (sojourn s)
  if (length(s_23)) {
    pD <- ggplot() +
      geom_histogram(aes(x = s_23, y = after_stat(density)),
                     binwidth = fd_binwidth(s_23),
                     color = "grey40", fill = "grey80") +
      geom_line(data = df_D, aes(s, f, color = "mixture f_S(s)"), linewidth = 1) +
      labs(title = "Ill → death elapsed time (sojourn s)", x = "s", y = "Density", color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  } else {
    pD <- ggplot() +
      xlim(0, 1) + ylim(0, 1) + theme_void(base_size = 11) +
      geom_text(aes(0.5, 0.5, label = "No ill→death times observed"))
  }
  
  # (E) Death time | path 1→2→3 (calendar time)
  if (length(death_123)) {
    pE <- ggplot() +
      geom_histogram(aes(x = death_123, y = after_stat(density)),
                     binwidth = fd_binwidth(death_123),
                     color = "grey40", fill = "grey80") +
      geom_line(data = df_E, aes(t, f, color = "theory f_T(t | 1→2→3)"), linewidth = 1) +
      labs(title = "Death time | path 1→2→3 (calendar time)", x = "t", y = "Density", color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  } else {
    pE <- ggplot() +
      xlim(0, 1) + ylim(0, 1) + theme_void(base_size = 11) +
      geom_text(aes(0.5, 0.5, label = "No 1→2→3 paths observed"))
  }
  
  # (F) Summary panel as a geom_text plot
  emp_P12 <- mean(got_ill)
  emp_P13 <- mean(dir_death)
  txtF <- sprintf(
    "Estimated path probabilities from simulation (n=%d):\n  P(1→2→3)  empirical = %.4f ; theory = %.4f\n  P(1→3)    empirical = %.4f ; theory = %.4f\n  Check: P12+P13 empirical = %.4f ; theory = %.4f",
    n, emp_P12, P12, emp_P13, P13, emp_P12 + emp_P13, P12 + P13
  )
  
  pF <- ggplot() +
    xlim(0, 1) + ylim(0, 1) + theme_void(base_size = 11) +
    geom_text(aes(x = 0, y = 1, label = txtF), hjust = 0, vjust = 1) +
    labs(title = "Summary: P12 & P13 shown below")
  
  # Compose 2×3 layout with patchwork
  my_plot <- (pA | pB | pC) /
    (pD | pE | pF)
  
  
  if (print_plot) {
    print(my_plot)
  }
  
  
  hazards = list(a12 = a12, a13 = a13, a23 = a23)
  class(hazards) <- c(class(hazards), "full_hazard")
  
  # 6) Return evaluators
  list(
    sim = sim,
    grid = list(t = t_grid, s = s_grid),
    # First-event
    pdf_first = stats::approxfun(t_grid, f_first, rule = 2),
    cdf_first = stats::approxfun(t_grid, F_first, rule = 2),
    # Illness time | illness
    pdf_ill_cond = stats::approxfun(t_grid, f12_cond, rule = 2),
    cdf_ill_cond = stats::approxfun(t_grid, F12_cond, rule = 2),
    # Direct death | direct
    pdf_dirdeath_cond = stats::approxfun(t_grid, f13_cond, rule = 2),
    cdf_dirdeath_cond = stats::approxfun(t_grid, F13_cond, rule = 2),
    # Sojourn s density and CDF
    pdf_ill2death_sojourn = stats::approxfun(s_grid, f23_mix, rule = 2),
    cdf_ill2death_sojourn = stats::approxfun(s_grid, F23_mix, rule = 2),
    # Death time | path 1->2->3 (calendar time)
    pdf_death_given_123 = stats::approxfun(t_grid, f_death_123, rule = 2),
    cdf_death_given_123 = stats::approxfun(t_grid, F_death_123, rule = 2),
    hazards = hazards,
    theory = list(
      t_grid = t_grid, S1 = S1, f12 = f12, f13 = f13, f_first = f_first,
      F12 = F12, F13 = F13, F_first = F_first, P12 = P12, P13 = P13,
      s_grid = s_grid, S23_sojourn = S23_mix
    ),
    plot = my_plot
  )
}



add_interval_censoring_to_illness <- function(dt, obs_interval = 1, obs_time_sd = 0.1) {
  # Required columns
  stopifnot(all(c("entry2", "death_time") %in% names(dt)))
  
  # Extract core times
  time_to_illness <- as.numeric(dt$entry2)
  time_to_death   <- as.numeric(dt$death_time)
  n <- length(time_to_illness)
  
  # ---- Censoring: subject-specific Uniform(0, 3 * death_i) with Inf-safe fallback
  finite_deaths_max <- max(time_to_death[is.finite(time_to_death)], 1)
  time_to_censor <- 1 + runif(n) * 3 * ifelse(is.finite(time_to_death), time_to_death, finite_deaths_max)
  
  # Observation cutoff per subject
  T_cutoff <- pmin(time_to_death, time_to_censor)
  
  # ---- Build an observation schedule up to the global max cutoff
  max_follow_up <- max(T_cutoff[is.finite(T_cutoff)], 0)
  grid <- seq(0, max_follow_up, by = obs_interval)
  
  obs_schedule <- matrix(rep(grid, n), nrow = n, byrow = TRUE)
  n_obs <- ncol(obs_schedule)
  
  # Jitter the schedule but preserve monotonicity
  if (n_obs > 1 && obs_time_sd > 0) {
    noise <- matrix(rnorm(n * n_obs, mean = 0, sd = obs_time_sd), nrow = n)
    clip <- max(min(obs_interval / 2 - 1e-6, 3 * obs_time_sd), 0)  # keep order
    if (clip > 0) {
      noise <- pmin(pmax(noise, -clip), clip)
      obs_schedule <- pmax(obs_schedule + noise, 0)
      # enforce increasing times row-wise
      obs_schedule <- t(apply(obs_schedule, 1, sort))
    }
  }
  # Force first visit at time 0
  obs_schedule[, 1] <- abs(obs_schedule[, 1])
  V_0 <- rep(0, n)
  
  # ---- Determine last healthy visit before illness (or before cutoff if earlier)
  # We use T_end = min(illness time, cutoff) row-wise
  T_end <- pmin(time_to_illness, T_cutoff)
  # Count visits strictly before T_end, row-wise
  idx_healthy <- rowSums(sweep(obs_schedule, 1, T_end, "<"))
  
  # V_healthy: if no visit before T_end, use baseline 0
  V_healthy <- V_0
  sel <- idx_healthy > 0
  V_healthy[sel] <- obs_schedule[cbind(which(sel), idx_healthy[sel])]
  
  # ---- First visit after the last healthy visit (candidate "ill" visit)
  idx_ill <- pmin(idx_healthy + 1L, n_obs)
  V_ill_candidate <- obs_schedule[cbind(seq_len(n), idx_ill)]
  
  # Illness occurs before cutoff?
  ill_before_cutoff <- is.finite(time_to_illness) & (time_to_illness <= T_cutoff)
  
  # Candidate ill-visit is valid only if it happens on/before cutoff and there is a next visit
  has_next_visit_before_cutoff <- ill_before_cutoff &
    (idx_healthy < n_obs) &
    (V_ill_candidate <= T_cutoff + 1e-12)
  
  V_ill <- ifelse(has_next_visit_before_cutoff, V_ill_candidate, NA_real_)
  
  # ---- Status classification (purely time-based; no reliance on 'path' labels)
  # 1: no illness observed (no V_ill), alive at cutoff
  # 2: no illness observed (no V_ill), died at cutoff
  # 3: interval-censored illness (have V_ill), alive at cutoff
  # 4: interval-censored illness (have V_ill), died at cutoff
  died_at_cutoff <- is.finite(time_to_death) & (time_to_death <= time_to_censor)
  has_interval   <- !is.na(V_ill)
  
  status <- integer(n)
  status[!has_interval & !died_at_cutoff] <- 1
  status[!has_interval &  died_at_cutoff] <- 2
  status[ has_interval & !died_at_cutoff] <- 3
  status[ has_interval &  died_at_cutoff] <- 4
  
  # For status 1 (alive, never observed ill), the illness time is right-censored at the last healthy visit
  T_obs <- T_cutoff
#  T_obs[status == 1] <- V_healthy[status == 1]
  
  # Return
  obs_dt <- data.table::data.table(
    V_0 = V_0,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status = factor(
      status,
      levels = 1:4,
      labels = c("healthy@cutoff (cens)", "died@cutoff (no illness observed)",
                 "interval illness, alive@cutoff", "interval illness, died@cutoff")
    )
  )
  
  # Store the (possibly wide) schedule as a list-column to avoid unintended column expansion
  schedule_list <- split(obs_schedule, row(obs_schedule))
  true_dt <- data.table::data.table(
    obs_schedule = schedule_list,
    time_to_censor = time_to_censor,
    time_to_illness = time_to_illness,
    time_to_death = time_to_death,
    T_cutoff = T_cutoff,
    status = obs_dt$status
  )
  
  list(obs = obs_dt, true = true_dt)
}

add_interval_censoring_to_illness_extension <- function(
    dt,
    scenario = 1L,                # 1..4: visit grids & missingness per Frydman–Szarek
    study_end = 1460,             # days (≈ 4 years)
    obs_time_sd = 100,            # days; N(0, sd) jitter at each planned visit (except baseline)
    seed = NULL
) {
  # Required columns
  stopifnot(all(c("entry2", "death_time") %in% names(dt)))
  n <- nrow(dt)
  time_to_illness <- as.numeric(dt$entry2)
  time_to_death   <- as.numeric(dt$death_time)
  
  # ---- Scenario-specific visit schedule (in days)
  # Scen. 1 & 3: every 6 months; Scen. 2 & 4: 0, 2, 4, 6 months; then every 6m through 24m; then yearly
  six_months <- 182.5
  base_schedule <- switch(
    as.character(scenario),
    "1" = seq(0, study_end, by = six_months),
    "3" = seq(0, study_end, by = six_months),
    "2" = {
      v <- c(0, 60, 120, 180, 360, 540, 720, 1080, 1460)
      v[v <= study_end]
    },
    "4" = {
      v <- c(0, 60, 120, 180, 360, 540, 720, 1080, 1460)
      v[v <= study_end]
    },
    stop("scenario must be 1, 2, 3, or 4")
  )
  m <- length(base_schedule)
  if (m < 2) stop("Study end too short for any post-baseline visits.")
  
  # ---- Build per-subject schedule with jitter; keep monotone and within [0, study_end]
  if (!is.null(seed)) set.seed(seed)
  obs_schedule <- matrix(rep(base_schedule, each = n), nrow = n, ncol = m)
  if (obs_time_sd > 0) {
    noise <- matrix(rnorm(n * (m - 1), mean = 0, sd = obs_time_sd), nrow = n, ncol = (m - 1))
    obs_schedule[, 2:m] <- obs_schedule[, 2:m] + noise
    obs_schedule[, 2:m] <- pmax(pmin(obs_schedule[, 2:m], study_end), 0)  # clamp
    # enforce increasing times row-wise; baseline fixed at 0
    obs_schedule <- t(apply(obs_schedule, 1, function(x) { x <- sort(x); x[1] <- 0; x }))
  }
  
  # ---- Missingness of the intermediate status (illness) by assessment; baseline never missing
  # p1 = 0.5% in Scen. 1–2; p1 = 3% in Scen. 3–4; doubles each assessment; cap at 0.98.
  p1 <- if (scenario %in% c(1L, 2L)) 0.005 else 0.03
  miss_prob <- pmin(0.98, p1 * 2^(0:(m - 2)))     # length m-1 (post-baseline)
  miss_mat  <- matrix(runif(n * (m - 1)) < rep(miss_prob, each = n), nrow = n, ncol = (m - 1))
  missing_status <- cbind(rep(FALSE, n), miss_mat)  # include baseline (never missing)
  
  # ---- Administrative end & death
  T_end <- pmin(time_to_death, study_end)
  died_by_end <- is.finite(time_to_death) & (time_to_death <= study_end)
  
  # ---- Determine interval [V_healthy, V_ill] if illness occurs before end and before death
  V_0 <- rep(0, n)
  V_healthy <- rep(0, n)
  V_ill <- rep(NA_real_, n)
  
  ill_before_end <- is.finite(time_to_illness) & (time_to_illness <= study_end) & (time_to_illness < time_to_death)
  
  for (i in seq_len(n)) {
    t_vis  <- obs_schedule[i, ]
    known  <- (!missing_status[i, ]) & (t_vis <= (T_end[i] + 1e-12))
    known[1] <- TRUE  # baseline known
    
    if (ill_before_end[i]) {
      idx_before <- which((t_vis < time_to_illness[i]) & known)
      j_star <- if (length(idx_before)) max(idx_before) else 1L
      idx_after <- which(seq_len(m) > j_star & known)
      if (length(idx_after)) {
        k <- min(idx_after)
        V_healthy[i] <- t_vis[j_star]
        V_ill[i]     <- t_vis[k]
      } else {
        # No known post-illness visit before T_end
        V_healthy[i] <- t_vis[j_star]
        V_ill[i]     <- NA_real_
      }
    } else {
      # No illness before end/death: last known healthy visit before T_end
      idx_known <- which(known & (t_vis <= (T_end[i] + 1e-12)))
      V_healthy[i] <- if (length(idx_known)) max(t_vis[idx_known]) else 0
      V_ill[i]     <- NA_real_
    }
  }
  
  # ---- Status categories (as in your original)
  status <- integer(n)
  status[is.na(V_ill) & !died_by_end] <- 1L
  status[is.na(V_ill) &  died_by_end] <- 2L
  status[!is.na(V_ill) & !died_by_end] <- 3L
  status[!is.na(V_ill) &  died_by_end] <- 4L
  
  # ---- Outputs
  obs_dt <- data.table::data.table(
    V_0 = V_0,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_end,
    status = factor(
      status, levels = 1:4,
      labels = c("healthy@admin_end (cens)",
                 "died@admin_end (no illness observed)",
                 "interval illness, alive@admin_end",
                 "interval illness, died@admin_end")
    )
  )
  
  # keep schedules and missingness as list-cols
  schedule_list <- split(obs_schedule, row(obs_schedule))
  missing_list  <- split(missing_status, row(missing_status))
  
  true_dt <- data.table::data.table(
    scenario = scenario,
    obs_schedule = schedule_list,
    missing_status = missing_list,
    time_to_illness = time_to_illness,
    time_to_death = time_to_death,
    study_end = study_end,
    died_by_end = died_by_end,
    status = obs_dt$status
  )
  
  structure(
    list(obs = obs_dt, true = true_dt),
    scenario = scenario,
    missing_p1 = p1,
    obs_time_sd = obs_time_sd,
    class = "fs_interval_obs"
  )
}

simulate_idm <- function(n, a12, a13, a23, print_plot = F) {
  res <- verify_illness_death(
    n,
    a12 = a12,
    a13 = a13,
    a23 = a23,
    print_plot = print_plot
  )
  
  res_w_ic <- add_interval_censoring_to_illness_extension(res$sim)
  list(
    obs = res_w_ic$obs,
    true_data_generation = res,
    true_censoring =  res_w_ic$true
  )
}


simulate_idm_constant_hazards <- function(
    n = 1000,
    a12 = 0.1,
    a13 = 0.2,
    a23 = 0.3, 
    print_plot = F) {
  
  a12_const <- function(t) rep(a12, length(t))
  a13_const <- function(t) rep(a13, length(t))
  a23_const <- function(s) rep(a23, length(s))
  simulate_idm(n, a12_const, a13_const, a23_const, print_plot = print_plot)
}
# res_constant <- simulate_idm_constant_hazards(n = 1000)

simulate_idm_weibull <- function(
    n = 1000,
    shape12 = 3, scale12 = 1,
    shape13 = 5, scale13 = 1,
    shape23 = 2, scale23 = 1,
    print_plot = F) {
  # validate inputs
  params <- c(shape12, scale12, shape13, scale13, shape23, scale23)
  if (any(!is.finite(params)) || any(params <= 0))
    stop("All shapes and scales must be positive and finite.")
  
  # Weibull hazard: (k/λ) * (t/λ)^(k-1), vectorized and safe at t=0
  h_weibull <- function(t, shape, scale) {
    t <- pmax(as.numeric(t), .Machine$double.eps)
    (shape / scale) * (t / scale)^(shape - 1)
  }
  
  a12 <- function(t) h_weibull(t, shape12, scale12)  # 1 -> 2
  a13 <- function(t) h_weibull(t, shape13, scale13)  # 1 -> 3
  a23 <- function(t) h_weibull(t, shape23, scale23)  # 2 -> 3 
  
  simulate_idm(n, a12, a13, a23, print_plot = print_plot)
}
# res_weibull <- simulate_idm_weibull()


simulate_idm_poly <- function(n = 1000, print_plot = F) {
  a12 <- function(x) dnorm(x)  # 1 -> 2
  a13 <- function(x) pmax(0.1,-0.12*(x-3)^4+1.2*(x-3))  # 1 -> 3
  a23 <- function(x) pmin(x^2,5)  # 2 -> 3 
  
  simulate_idm(n, a12, a13, a23, print_plot = print_plot)
}
# res_poly <- simulate_idm_poly()


