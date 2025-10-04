#### Interval manipulation ####
##### intersect ####
# generic
intersect <- function(x, y, ...) {
  if (inherits(x, "interval")) {
    UseMethod("intersect", x)
  } else if (inherits(y, "interval")) {
    UseMethod("intersect", y)
  } else {
    base::intersect(x, y, ...)
  }
}

intersect.default <- function(x, y, ...) base::intersect(x, y, ...)

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

##### get_interval ####

# generic
get_interval <- function(x, ...) {
  UseMethod("get_interval")
}

get_interval.matrix <- function(m, L_open = F, R_open = F) {
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


##### as.interval ####

as.interval <- function(x, L_open = F, R_open = F) {
  if (ncol(x) != 2) stop("LR mat must be of dim m x 2")
  type = paste(ifelse(L_open, "o", "c"), ifelse(R_open, "o", "c"), sep = "_")
  
  class(x) <- c("interval", type, class(x))
  x
}

##### contains ####

# generic
is_subset <- function(A, B, ...) {
  UseMethod("is_subset")
}

# checks if A ⊂ B
is_subset.interval <- function(A, B, strict_subset = F) {
  
  if(!(inherits(B, "interval") & all(dim(A) == dim(B)))) stop("B should be an of same dim as A interval")

  if (inherits(B, "c_c")) {
    l_compare <- `<=`
    r_compare <- `<=`
  } else if (inherits(B, "c_o")) {
    l_compare <- `<=`
    r_compare <- `<`
  } else if (inherits(B, "o_c")) {
    l_compare <- `<`
    r_compare <- `<=`
  } else if (inherits(B, "o_o")) {
    l_compare <- `<`
    r_compare <- `<`
  }
  
  if (inherits(A, "c_c")) {
  } else if (inherits(A, "c_o")) {
    r_compare <- `<=`
  } else if (inherits(A, "o_c")) {
    l_compare <- `<=`
  } else if (inherits(A, "o_o")) {
    l_compare <- `<=`
    r_compare <- `<=`
  }
  
  return(l_compare(B[,1],A[,1]) & r_compare(A[,2],B[,2]))
}




#### From paper specific ####

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

# Comment: the intervals input determines if the interval is open or closed
product_over_t_stars <- function(intervals, T_star, lambda_n) {
  if(!inherits(intervals, "interval")) stop("intervals need to be of type interval")
  T_stars_to_prod_over <- intersect(T_star, intervals)
  prod_lambdas <- lambda_n[which(T_star %in% T_stars_to_prod_over)]
  prod(1-prod_lambdas)
}

product_over_t_stars_one_interval <- function(L, R, L_open, R_open,T_star, lambda_n) {
  intervals <- as.interval(matrix(c(L, R), ncol = 2), L_open, R_open)
  T_stars_to_prod_over <- intersect(T_star, intervals)
  prod_lambdas <- lambda_n[which(T_star %in% T_stars_to_prod_over)]
  prod(1-prod_lambdas)
}



