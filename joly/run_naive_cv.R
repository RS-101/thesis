# USING INTERCEPT = F IN SPLINES. IF WE WANT INTERCEPT WE NEED TO CHANGE SPLINE DESIGN FOR PEN MAT


library(tidyverse)

debugSource("likelihood_hazard_splines.R")


## ---------------------------
## Helper: fit with penalty
## ---------------------------
do_likelihood_optim <- function(sim_dat, n_knots, degree = 3, penalizer, knots = NULL) {
  n_par <- n_knots + degree - 1
  
  # Use training-derived knots (important for CV fairness). Allow override.
  if (is.null(knots)) {
    knots_a01 <- seq(min(sim_dat$V_0),  max(sim_dat$T_obs), length.out = n_knots)
    knots_a02 <- knots_a01
    knots_a12 <- knots_a01
    knots <- list(a01 = knots_a01, a02 = knots_a02, a12 = knots_a12)
  } else {
    stopifnot(length(knots$a01) == n_knots, length(knots$a02) == n_knots, length(knots$a12) == n_knots)
  }
  
  # penalty matrix for each hazard (must exist in your env)
  P01 <- pen_mat_m_splines(knots$a01)
  P02 <- pen_mat_m_splines(knots$a02)
  P12 <- pen_mat_m_splines(knots$a12)
  
  # NOTE: We MINIMIZE:  -loglik + lambda * penalty_sum
  # Your original code subtracted the penalty and also used '-' between hazards,
  # which effectively rewards roughness. We fix both issues here.
  obj_fun <- function(x) {
    theta <- list(
      a01 = x[1:n_par],
      a02 = x[(n_par+1):(2*n_par)],
      a12 = x[(2*n_par+1):(3*n_par)]
    )
    ll <- full_log_likehood(
      V_0      = sim_dat$V_0,
      V_healthy= sim_dat$V_healthy,
      V_ill    = sim_dat$V_ill,
      T_obs    = sim_dat$T_obs,
      status   = sim_dat$status,
      theta    = theta,
      degree   = degree,
      knots    = knots
    )$loglik
    
    pen <- drop(t(theta$a01) %*% P01 %*% theta$a01) +
      drop(t(theta$a02) %*% P02 %*% theta$a02) +
      drop(t(theta$a12) %*% P12 %*% theta$a12)
    
    val <- -(ll - penalizer * pen)
    if (!is.finite(val)) return(1e10)
    val
  }
  
  x0 <- rep(1, 3*n_par)
  # Warm check
  invisible(obj_fun(x0))
  
  res <- optim(
    par   = x0,
    fn    = obj_fun,
    method= "L-BFGS-B",
    lower = rep(0.0001, 3*n_par),
    upper = rep(13,     3*n_par)
  )
  
  theta_hat <- list(
    a01 = res$par[1:n_par],
    a02 = res$par[(n_par+1):(2*n_par)],
    a12 = res$par[(2*n_par+1):(3*n_par)]
  )
  
  fit_train <- full_log_likehood(
    V_0      = sim_dat$V_0,
    V_healthy= sim_dat$V_healthy,
    V_ill    = sim_dat$V_ill,
    T_obs    = sim_dat$T_obs,
    status   = sim_dat$status,
    theta    = theta_hat,
    degree   = degree,
    knots    = knots
  )
  
  # Return everything needed for validation and later use
  list(
    theta        = theta_hat,
    knots        = knots,
    degree       = degree,
    train_loglik = fit_train$loglik,
    train_pen    = {
      P01 <- pen_mat_m_splines(knots$a01)
      P02 <- pen_mat_m_splines(knots$a02)
      P12 <- pen_mat_m_splines(knots$a12)
      drop(t(theta_hat$a01) %*% P01 %*% theta_hat$a01) +
        drop(t(theta_hat$a02) %*% P02 %*% theta_hat$a02) +
        drop(t(theta_hat$a12) %*% P12 %*% theta_hat$a12)
    },
    train_obj    = -fit_train$loglik + penalizer * (
      drop(t(theta_hat$a01) %*% P01 %*% theta_hat$a01) +
        drop(t(theta_hat$a02) %*% P02 %*% theta_hat$a02) +
        drop(t(theta_hat$a12) %*% P12 %*% theta_hat$a12)
    ),
    optim        = res,
    hazards      = fit_train$hazards # passthrough for convenience
  )
}

## ---------------------------
## Helper: evaluate LL on new data
## ---------------------------
loglik_on_data <- function(sim_dat, theta, degree, knots) {
  out <- full_log_likehood(
    V_0      = sim_dat$V_0,
    V_healthy= sim_dat$V_healthy,
    V_ill    = sim_dat$V_ill,
    T_obs    = sim_dat$T_obs,
    status   = sim_dat$status,
    theta    = theta,
    degree   = degree,
    knots    = knots
  )
  out$loglik
}

## ---------------------------
## K-fold CV across (n_knots, penalizer) grid
## ---------------------------
cv_penalized_splines <- function(sim_dat,
                                 n_knots_grid,
                                 penalizer_grid,
                                 K = 5,
                                 degree = 3,
                                 seed = 1,
                                 verbose = TRUE) {
  stopifnot(K >= 2)
  set.seed(seed)
  n <- nrow(sim_dat)
  folds <- sample(rep(1:K, length.out = n))
  
  res_rows <- length(n_knots_grid) * length(penalizer_grid)
  results <- vector("list", res_rows)
  row_id <- 0L
  
  for (nk in n_knots_grid) {
    for (lam in penalizer_grid) {
      row_id <- row_id + 1L
      if (verbose) message(sprintf("[CV] n_knots=%d, lambda=%.6g", nk, lam))
      
      fold_ll  <- numeric(K)
      fold_n   <- numeric(K)
      train_ll <- numeric(K)
      
      for (k in 1:K) {
        idx_val <- which(folds == k)
        idx_tr  <- setdiff(seq_len(n), idx_val)
        
        dat_tr  <- sim_dat[idx_tr, , drop = FALSE]
        dat_val <- sim_dat[idx_val, , drop = FALSE]
        
        # Fit on training fold -> returns theta and training-derived knots
        fit_k <- do_likelihood_optim(dat_tr, n_knots = nk, degree = degree, penalizer = lam)
        
        # Evaluate out-of-fold loglik using training knots
        ll_val <- try(loglik_on_data(dat_val, theta = fit_k$theta, degree = degree, knots = fit_k$knots),
                      silent = TRUE)
        if (inherits(ll_val, "try-error") || !is.finite(ll_val)) ll_val <- -Inf
        
        fold_ll[k]  <- ll_val
        fold_n[k]   <- nrow(dat_val)
        train_ll[k] <- fit_k$train_loglik
      }
      
      # Aggregate: total OOF loglik and per-subject mean for comparability
      total_ll_oof <- sum(fold_ll)
      mean_ll_oof  <- sum(fold_ll) / n
      mean_ll_tr   <- sum(train_ll) / n
      
      results[[row_id]] <- data.frame(
        n_knots          = nk,
        lambda           = lam,
        oof_loglik_total = total_ll_oof,
        oof_loglik_mean  = mean_ll_oof,
        train_loglik_mean= mean_ll_tr,
        K                = K,
        stringsAsFactors = FALSE
      )
    }
  }
  
  cv_table <- do.call(rbind, results)
  # pick the model that maximizes out-of-fold log-likelihood
  best_idx <- which.max(cv_table$oof_loglik_mean)
  best <- cv_table[best_idx, , drop = FALSE]
  
  # Refit on full data with best hyperparameters
  final_fit <- do_likelihood_optim(sim_dat,
                                   n_knots  = best$n_knots,
                                   degree   = degree,
                                   penalizer= best$lambda)
  
  list(
    cv_table  = cv_table[order(-cv_table$oof_loglik_mean), ],
    best      = best,
    final_fit = final_fit,
    folds     = folds
  )
}

## ---------------------------
## Example usage
## ---------------------------
# Choose your grids (adjust as needed)
# n_knots must be >= 3; include a range to let CV choose smoothness via basis size + lambda
# penalizer grid on log-scale is typical
grid_nk <- c(3, 4, 5, 7, 9)
grid_lam <- 10^seq(-4, 4, length.out = 9)
set.seed(42)
cv_res <- cv_penalized_splines(res_w_ic$obs, n_knots_grid = grid_nk, penalizer_grid = grid_lam, K = 5)
cv_res$cv_table
cv_res$final_fit  # contains theta, knots, hazards, etc.

res_full <- cv_res$final_fit

# Collect your functions
funs <- list(
  a01 = res_full$hazards$a01,
  a02 = res_full$hazards$a02,
  a12 = res_full$hazards$a12,
  true_a01 = res$hazards$a12,
  true_a02 = res$hazards$a13,
  true_a12 = res$hazards$a23
)

# Grid of x values
xvals <- seq(0, 13, length.out = 200)



# Evaluate each function
df <- lapply(names(funs), function(nm) {
  data.frame(
    x = xvals,
    y = sapply(xvals, funs[[nm]]),
    fun = nm
  )
}) %>% bind_rows()

df <- df %>% mutate(type = ifelse(str_detect(fun, "true"), "true", "estimate"))


# Plot
ggplot(df, aes(x = x, y = y, colour = fun)) +
  geom_line(size = 1, aes(linetype = type)) +
  labs(title = "Spline hazard",
       x = "x", y = "value") +
  theme_minimal()


