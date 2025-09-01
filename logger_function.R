create_logger <- function() {
  log_env <- new.env()
  log_env$records <- list()
  
  log_env$log <- function(iter, loss, metric, extra = list()) {
    log_env$records[[length(log_env$records) + 1]] <- c(
      list(
        iteration = iter,
        loss = loss
      ),
      extra
    )
  }
  
  log_env$get_log <- function() {
    dplyr::bind_rows(log_env$records)
  }
  
  return(log_env)
}
