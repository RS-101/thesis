library(data.table)

simulate_time_to <- function(n = 1000, censoring_time = 100) {
  # parameters
  mean_time_01 <- 3 # Ill from state 0
  mean_time_02 <- 10 # Death from state 0
  mean_time_12 <- 5
  
  a01 <- function(x) 1/mean_time_01
  a02 <- function(x) 1/mean_time_02
  a12 <- function(x) 1/mean_time_12
  
  # latent event times
  T01 <- rexp(n, 1/mean_time_01)
  T02 <- rexp(n, 1/mean_time_02)
  T12 <- rexp(n, 1/mean_time_12)
  
  path <- rep(NA_character_, n)
  
  path[pmin(T01, T02) >= censoring_time] <- "00"
  
  # die before illness
  path[T02 < T01 & T02 < censoring_time] <- "02"
  
  # ill before censoring; then either die before censoring or be censored in 1
  ill <- (T01 <= T02) & (T01 < censoring_time)
  path[ill & (T01 + T12 < censoring_time)]  <- "012"
  path[ill & (T01 + T12 >= censoring_time)] <- "01"

  time_to_illness = ifelse(path == "01" | path == "012", T01, NA)
  time_to_death = ifelse(path == "012", T01+T12, T02)
  time_to_death[time_to_death > censoring_time] = NA
  time_to_censor = ifelse(path == "00" | path == "01", censoring_time, NA) 

  list(
    data = data.table(
      time_to_illness = time_to_illness,
      time_to_death = time_to_death,
      time_to_censor = time_to_censor,
      path
    ), 
    hazards = list(
      "a01" = a01,
      "a02" = a02,
      "a12" = a12
    )
  )
}





add_interval_censoring_to_illness <- function(dt, obs_interval = 1, obs_time_sd = 0.1) {
  
  
  time_to_illness <- dt$time_to_illness
  time_to_censor <- dt$time_to_censor
  time_to_death <- dt$time_to_death
  path <- dt$path
  n <- length(time_to_illness)
  
  max_follow_up <- max(time_to_censor, time_to_death, na.rm = T)
  obs_schedule <- matrix(
    rep(seq(0, max_follow_up, by = obs_interval), n), 
    nrow = n, byrow = TRUE
  )
  
  n_obs <- ncol(obs_schedule)
  
  # Truncated normal noise distribution
  obs_noise <- matrix(
    rnorm(n * n_obs, mean = 0, sd = obs_time_sd),
    nrow = n
  )
  obs_noise <- pmax(pmin(obs_noise, (obs_interval - 0.1)), -(obs_interval + 0.1))
  obs_schedule <- pmax(obs_schedule + obs_noise, 0)
  
  V_0 <- obs_schedule[, 1] <- rep(0, n)
  
  idx_healthy <- ifelse(!(path == "00" | path == "02"), rowSums(obs_schedule < time_to_illness), rowSums(obs_schedule < max_follow_up))
  
  V_healthy <- obs_schedule[cbind((1:n), idx_healthy)]
  
  V_ill <- ifelse(path == "01" | path == "012",
                  obs_schedule[cbind(1:n, pmin(idx_healthy + 1, n_obs))],
                  NA)
  
  V_m <- pmax(V_ill, V_healthy, na.rm = T)
  
  T_obs <- pmin(time_to_censor, time_to_death, na.rm = T)
  
  # Subject is healthy at Vm and still alive at T (no information about illness at T )
  status <- ifelse((path == "00" | path == "01") & V_m == V_healthy, 1, NA)
  
  # Subject is healthy at Vm and die at T (we do not know if subject is ill or healthy at the time of subjects death)
  status <- ifelse((path == "02" | path == "012") & V_m == V_healthy, 2, status)
  
  # Subject is healthy at Vk, (k < m), ill at Vk+1, and still alive at T
  status <- ifelse(path == "01" & V_healthy < V_m , 3, status)
  
  # Subject is healthy at Vk, (k < m), ill at Vk+1, and dies at T :
  status <- ifelse(path == "012" & V_healthy < V_m, 4, status)
  
  T_obs[status == 1] = V_healthy[status == 1]
  
  list(
    obs = data.table(
      V_0 = V_0, 
      V_healthy = V_healthy, 
      V_ill = V_ill,
      T_obs = T_obs, 
      status = factor(status)
      ),
    true = data.table(
      obs_schedule = obs_schedule, 
      time_to_censor = time_to_censor,
      time_to_illness = time_to_illness,
      time_to_death = time_to_death, 
      status = factor(status)
      )
    )
}


simulate_data <- function(n = 1000, censoring_time = 1) {
  res <- simulate_time_to(n, censoring_time)
  # print(res)
  list(data = add_interval_censoring_to_illness(res$data), 
       hazards = res$hazards)
}


# 
# my_data <- simulate_data(1000, censoring_time = 3)
# my_data$obs$status %>% as.factor() %>% summary()
