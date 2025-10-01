source("joly/functions_likelihood.R")

#### Perform K-fold CV ####
cv_penalized_splines <- function(sim_dat,
                                 n_knots_grid,
                                 penalizer_a12,
                                 penalizer_a13,
                                 penalizer_a23,
                                 K = 5,
                                 degree = 3,
                                 seed = 1,
                                 verbose = TRUE) {
  stopifnot(K >= 2)
  set.seed(seed)
  n <- nrow(sim_dat)
  folds <- sample(rep(1:K, length.out = n))
  
  
  pen_vals <- expand.grid(penalizer_a12, penalizer_a13, penalizer_a23)
  rows_pen_vals <- nrow(pen_vals)
  
  res_rows <- length(n_knots_grid) * rows_pen_vals
  results <- vector("list", res_rows)
  row_id <- 0L
  
  for (nk in n_knots_grid) {
    for (lam in seq_len(rows_pen_vals)) {
      row_id <- row_id + 1L
      if (verbose) {
        message(sprintf(
          "[CV] n_knots=%d, lambda=(%.6g, %.6g, %.6g)",
          nk,
          pen_vals[lam, 1],
          pen_vals[lam, 2],
          pen_vals[lam, 3]
        ))
      }
      
      fold_ll  <- numeric(K)
      fold_n   <- numeric(K)
      train_ll <- numeric(K)
      
      for (k in 1:K) {
        message(sprintf(
          "[CV] fold %d of %.6g",
          k,K
        ))
        idx_val <- which(folds == k)
        idx_tr  <- setdiff(seq_len(n), idx_val)
        
        dat_tr  <- sim_dat[idx_tr, , drop = FALSE]
        dat_val <- sim_dat[idx_val, , drop = FALSE]
        
        # Fit on training fold -> returns theta and training-derived knots
        fit_k <- do_likelihood_optim(dat_tr,
                                     n_knots = nk,
                                     degree = degree,
                                     penalizer = as.numeric(pen_vals[lam,]))
        
        # Evaluate out-of-fold loglik using training knots
        fold_ll_val <- full_log_likehood(
          V_0       = dat_val$V_0,
          V_healthy = dat_val$V_healthy,
          V_ill     = dat_val$V_ill,
          T_obs     = dat_val$T_obs,
          status    = dat_val$status,
          theta     = fit_k$theta,
          degree    = degree,
          knots     = fit_k$knots
        )$loglik
        

        fold_ll[k]  <- fold_ll_val
        fold_n[k]   <- nrow(dat_val)
        train_ll[k] <- fit_k$loglik
      }
      
      # Aggregate: total OOF loglik and per-subject mean for comparability
      total_ll_oof <- sum(fold_ll)
      mean_ll_oof  <- total_ll_oof / n
      mean_ll_tr   <- sum(train_ll) / n
      
      results[[row_id]] <- data.frame(
        n_knots          = nk,
        lambda_12           = pen_vals[lam, 1],
        lambda_13           = pen_vals[lam, 2],
        lambda_23           = pen_vals[lam, 3],
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
                                   penalizer= c(best$lambda_12,
                                                best$lambda_13,
                                                best$lambda_23) )
  
  list(
    cv_table  = cv_table[order(-cv_table$oof_loglik_mean), ],
    best      = best,
    final_fit = final_fit,
    folds     = folds
  )
}
