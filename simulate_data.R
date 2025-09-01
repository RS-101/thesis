library(data.table)

simulate_time_to <- function(n = 1000, censoring_time = 100) {
  # parameters
  mean_time_01 <- 3 # Ill from state 0
  mean_time_02 <- 6 # Death from state 0
  mean_time_12 <- 1
  
  # latent event times
  T01 <- rexp(n, 1/mean_time_01)
  T02 <- rexp(n, 1/mean_time_02)
  T12 <- rexp(n, 1/mean_time_12)
  

  
  path <- ifelse(T01 > T02 & T02 < censoring_time,
                 "02",
                 ifelse(T01 > T02 & T02 >= censoring_time,
                        "00",
                         ifelse((T01 + T12) < censoring_time,
                                "012",
                                "01")))


  time_to_illness = ifelse(path == "01" | path == "012", T01, NA)
  time_to_death = ifelse(path == "012", T01+T12, T02)
  time_to_death[time_to_death > censoring_time] = NA
  time_to_censor = ifelse(path == "00", censoring_time, NA) 

  data.table(
    time_to_illness = time_to_illness,
    time_to_death = time_to_death,
    time_to_censor = time_to_censor,
    path
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
  idx_healthy <- ifelse(!(path == "00" | path == "02"), rowSums(obs_schedule < time_to_illness), n_obs)
  
  V_healthy <- obs_schedule[cbind((1:n), idx_healthy)]
  
  V_ill <- ifelse(path == "01" | path == "012",
                  obs_schedule[cbind(1:n, pmin(idx_healthy + 1, n_obs))],
                  NA)
  
  # idx_death <- ifelse(!(path == "00" | path == "01"),
  #                     rowSums(obs_schedule < time_to_death),
  #                     NA)
  # 
#  last_before_death <- obs_schedule[cbind(1:n, pmin(idx_death, n_obs))]
  
  T_obs <- pmin(time_to_censor, time_to_death, na.rm = T)
  
  # Subject is healthy at Vm and still alive at T (no information about illness at T )
  status <- ifelse(path == "00", 1, NA)
  # Subject is healthy at Vm and die at T (we do not know if subject is ill or healthy at the time of subjects death)
  status <- ifelse(path == "02" | (V_healthy <= T_obs & T_obs <= V_ill), 2, status)
  # Subject is healthy at Vk, (k < m), ill at Vk+1, and still alive at T
  status <- ifelse(path == "01" !is.na(V_ill) & (V_ill <= T_obs & time_to_death != T_obs), 3, status) 
  # Subject is healthy at Vk, (k < m), ill at Vk+1, and dies at T :
  status <- ifelse(!is.na(V_ill) & (V_ill < T_obs & T_obs == time_to_death), 4, status)
  
  V_healthy[status == 1] = T_obs[status == 1]
  V_ill[status == 1] = NA  
  
  V_ill[status == 2] = NA  
  
  list(
    obs = data.table(
      V_0 = V_0, 
      V_healthy = V_healthy, 
      V_ill = V_ill,
      T_obs = T_obs, 
      status = status
      ),
    true = data.table(
      obs_schedule = obs_schedule, 
      time_to_last_obs = time_to_censor,
      time_to_illness = time_to_illness,
      time_to_death = time_to_death, 
      status = status
      )
    )
}


simulate_data <- function(n = 1000, censoring_time = 1) {
  res <- simulate_time_to(n, censoring_time)
  print(res)
  add_interval_censoring_to_illness(res)
}



my_data <- simulate_data(100)
my_data$obs
