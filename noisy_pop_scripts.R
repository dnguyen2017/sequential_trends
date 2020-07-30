# simulation code for populations with process and observation error
library(ggplot2)

sim_pop <- function (xinit,     # initial log population size
                        lambda, # population growth rate
                        b,      # density dependence
                        phi,    # process noise autocorrelation [-1,1]
                        sd_proc, # process noise: N(0, sd_proc)
                        sd_obs,  # observation noise: N(0, sd_obs)
                        tfinal,
                        nsim) {
  
  # init storage list for all simulation
  sim_list <- vector("list", length = nsim)
  
  for (i in seq_along(sim_list)) {
    # init storage for current simulation
    log_x <- vector("numeric", length = tfinal)
    log_y <- vector("numeric", length = tfinal)
    
    proc_error <- rnorm(n = 1, mean = 0, sd = sd_proc)
    
    log_x[1] <- xinit  
    log_y[1] <- xinit + rnorm(n = 1, mean = 0, sd = sd_obs) 
    
    # sim population dynamics
    for (j in 2:tfinal) {
      proc_error <- phi * proc_error + rnorm(n = 1, mean = 0, sd = sd_proc)
      log_x[j] <- lambda + b * log_x[j-1] + proc_error
      log_y[j] <- log_x[j] + rnorm(n = 1, mean = 0, sd = sd_obs)
    }
    # save current simulation
    sim_list[[i]] <- (tidyr::tibble(time = 1:tfinal,
                                    log_x = log_x,
                                    log_y = log_y))
  } 
  
  return(dplyr::bind_rows(sim_list, .id = "sim"))
}

set_params_expand <- function (xinit = c(), lambda = c(), b = c(), phi = c(), sd_proc = c(), sd_obs = c()) {
  return(tidyr::expand_grid(xinit = xinit, lambda = lambda, b = b, phi = phi, sd_proc = sd_proc, sd_obs = sd_obs))
}

set_params <- function (xinit, lambda, b, phi, sd_proc, sd_obs) {
  return(tidyr::tibble(xinit = xinit, lambda = lambda, b = b, phi = phi, sd_proc = sd_proc, sd_obs = sd_obs))
}

multi_sim_pop <- function (params, tfinal, nsim) {
  sim_list <- vector("list", length = nrow(params))
  
  for ( i in 1:nrow(params)) {
    sim_list[[i]] <- with(params[i,],
                          sim_pop(xinit = xinit, lambda = lambda, b = b, phi = phi, sd_proc = sd_proc, sd_obs = sd_obs,
                             tfinal = tfinal, nsim = nsim)
                            )
    # add char col with param values
    par_label <- paste(names(params[i,]), "=", round(params[i,], 2), collapse = ",")
    sim_list[[i]]$par <- rep(par_label, times = nrow(sim_list[[i]]))  #rep(paste(round(params[i,], 2), collapse = ", "), times = nrow(sim_list[[i]]))
  }
  return(dplyr::bind_rows(sim_list))
}

       