# Clean version with only the estimator functions (no test code)

calc_F_12 <- function(z_i, Q_i, I) {
  l_1 = Q_i[1,1]
  r_I = Q_i[I,2]
  
  z_i <- z_i[1:I]
  
  raw_F_12 <- Vectorize(function(s) {
    if(s <= l_1) {
      return(0)
    } else if (s < r_I) {
      p <- which(Q_i[1:(I-1),2] <= s & s < Q_i[2:I,1])
      if (length(p) == 0) return(NA)
      
      return(sum(z_i[1:p]))
    } else {
      return(sum(z_i)) # sum up till I
    }
  }, vectorize.args = "s")
  
  fill_F_12 <- Vectorize(function(s) {
    if(s <= l_1) {
      return(0)
    } else if (s < r_I) {
      p <- min(which(s < Q_i[2:I,1]))
      return(sum(z_i[1:p]))
    } else {
      return(sum(z_i)) # sum up till I
    }
  }, vectorize.args = "s")
  
  list(
    raw = raw_F_12,
    fill = fill_F_12
  )
}

calc_F_13 <- function(z_i, Q_i_mark, I, I_mark) {
  
  z_i <- z_i[(I+1):I_mark]
  Vectorize(function(s) {
    sum(z_i * as.numeric(Q_i_mark <= s))
  }, vectorize.args = "s")
}

calc_F <- function(z_i, Q_i, Q_i_mark, I, I_mark) {
  F12 <- calc_F_12(z_i, Q_i, I)$fill
  F13 <- calc_F_13(z_i, Q_i_mark, I, I_mark)
  
  function(s) {
    F12(s)+F13(s)
  }
}

calc_A_12 <- function(z_i, Q_i, Q_i_mark, I, I_mark) {
  
  
  F_hat <- calc_F(z_i, Q_i, Q_i_mark, I, I_mark)
  
  F_hat_at_l <- F_hat((Q_i[,1]-1e-10))
  z_i <- z_i[1:I]
  
  Vectorize(function(s) {
    index <- (Q_i[,2] <= s)  
    sum(z_i[index]*1/(1-F_hat_at_l[index]))
  })
}

calc_A_13 <- function(z_i, Q_i, Q_i_mark, I, I_mark) {
  
  
  F_hat <- calc_F(z_i, Q_i, Q_i_mark, I, I_mark)
  
  F_hat_at_l <- F_hat((Q_i_mark-1e-10))
  z_i <- z_i[(I+1):I_mark]
  
  Vectorize(function(s) {
    index <- (Q_i_mark <= s)  
    sum(z_i[index]*1/(1-F_hat_at_l[index]))
  })
}

calc_A_23 <- function(lambda_n, t_star_n) {
  Vectorize(function(t) {
    index <- (t_star_n <= t)  
    sum(lambda_n[index])
  })
}
