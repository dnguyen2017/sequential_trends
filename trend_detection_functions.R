# linear trend detection

### simulation functions

# linear population trend model w/ normal process error
sim_lin_proc <- function (x0, t, trend, sd) {
  # allow for constant trend or vector of varying trends (automatically recycles vector if lenght(trend) < t)
  trend_val = rep(trend, length = t)
  sd_val = rep(sd, length = t)
  
  pop <- vector("numeric", length = t)
  pop[1] <- x0
  for (i in 1:(length(pop) - 1)) {
    pop[i + 1] <- rnorm(1, pop[i] + trend_val[i], sd_val[i])
  }
  return(tibble(time = 1:t, pop = pop, trend = trend_val, sd = sd_val))
}

# add normal and unbiased observation noise to a population time series
obs_model <- function (x, sd) {
  obs <- x + rnorm(n = length(x), mean = 0, sd = sd)
  return(obs)
}

# simulate sim_lin_proc nsim times using common parameters
multi_sim <- function(nsim, x0, t, trend, sd) {
  simulations <- lapply(1:nsim, function(x) sim_lin_proc(x0, t, trend, sd))
  return(bind_rows(simulations, .id = "simulation") %>% mutate(simulation = as.numeric(simulation)))
}

### analysis functions

# running calculation for mle estimates of normal R.V.s

run_stat <- function(x) {
  # add checks here
  
  # init vectors for statistics
  m_ <- vector("numeric", length = length(x)) 
  s_ <- s2_mle <- sd_mle <- s2_ <- sd_ <- m_
  
  # init value of mean is x1. All other init values are 0
  m_[1] <- x[1]
  
  # calculate running statistics
  for (i in 2:length(x)) {
    m_[i] <- m_[i-1] + (x[i] - m_[i-1])/i
    s_[i]  <- s_[i-1] + (x[i] - m_[i-1]) * (x[i] - m_[i])
    s2_mle[i] <- s_[i]/(i)
    sd_mle[i] <- sqrt(s2_mle[i])
    s2_[i] <- s_[i]/(i-1)
    sd_[i] <- sqrt(s2_[i])
  }
  
  # return: id = observation number, obs = observation value, mu_mle, sd_mle, sd_corrected
  return(tidyr::tibble(x, mu_mle = m_, sd_mle, sd_corrected = sd_)) # %>% tibble::rowid_to_column(var = "id"))
}

#out %>% ggplot() + geom_line(aes(x = id, y = mu_mle))
run_stat(rnorm(1000))
