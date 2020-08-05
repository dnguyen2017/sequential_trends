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

# apply FUNC cumulatively when NA are present
# https://stackoverflow.com/questions/25576358/calculate-cumsum-while-ignoring-na-values/25576972
cumSkipNA <- function(x, FUNC)
{
  d <- deparse(substitute(FUNC))
  funs <- c("max", "min", "prod", "sum")
  stopifnot(is.vector(x), is.numeric(x), d %in% funs)
  FUNC <- match.fun(paste0("cum", d))
  x[!is.na(x)] <- FUNC(x[!is.na(x)])
  x
}

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
  return(tidyr::tibble(diff_1 = x, mu_mle = m_, sd_mle, sd_corrected = sd_)) # %>% tibble::rowid_to_column(var = "id"))
  #return(tidyr::tibble(diff_1 = c(NA,x), mu_mle = c(NA,m_), sd_mle = c(NA,sd_mle), sd_corrected = c(NA,sd_) )) # %>% tibble::rowid_to_column(var = "id"))
}

#out %>% ggplot() + geom_line(aes(x = id, y = mu_mle))
#run_stat(rnorm(1000))

## helpers for seq_lm
get_pval <- function(model) {
  return(summary(model)$coefficients[2,4])
}
get_slope <- function(model) {
  return(model$coefficients[2])
}

# fixed sample size method for trend
seq_lm <- function (data, nsim) {
  # pre-allocate a list to store dfs
  lm_list <- vector("list", length = nsim)
  
  # fit lm and extract slope estimates and p-values
  for (i in seq_along(lm_list)) {
    # create empty cols for pval and trend estimate
    data$pval <- NA
    data$trend_est <- NA
    
    # grab current time series
    current_ts <- data %>% filter(simulation == i)
    
    # fit lm to all subsets t \in {{1,2}, {1,2,3}, {1,2,3,4} ...} in current_ts and extract slope est and pval
    for (j in 2:nrow(current_ts)) {
      lm_fit <- lm(pop ~ time, data = current_ts[1:j,])
      current_ts[j, "trend_est"] <- get_slope(lm_fit)
      current_ts[j, "pval"] <- get_pval(lm_fit)
    }
    # save calculated values
    lm_list[[i]] <- current_ts
    
  }
  
  return(bind_rows(lm_list))
}

# one sided sprt calculation
objfn <- function(trend, x1, x0, sd) {
  # likelihood( x1 | x0, r, dispersion )
  dnorm(x = x1, mean = x0 + trend, sd = sd)
}


#  find restricted MLE under h_0 or h_1
likelihood <- function (x_prev, x_now, sd, interval = c(), guess = 1.5) {
  opt <- optim(par = guess, fn = objfn, control = list(fnscale = -1), 
               method = "Brent", lower = interval[1], upper = interval[2],
               x0 = x_prev, x1 = x_now, sd = sd)
  return(opt$value)
}

# sequential test functions
calc_sr <- function(data, accept_reg = c(0,0), reject_reg = c(0, 5), sd_est) {
  
  # initialize empty columns for likelihood calculations
  data$lik0 <- NA
  data$lik1 <- NA
  
  # parameter space
  theta0 <- accept_reg
  theta1 <- reject_reg
  
  # use midpoint of interval as init guess for optim (Brent's method)
  guess0 <- accept_reg[1] + (accept_reg[2] - accept_reg[1]) / 2 
  guess1 <- reject_reg[1] + (reject_reg[2] - reject_reg[1]) / 2 
  
  sim_number <- 1
  for (i in 2:nrow(data)) {
    
    if (data[i, "simulation"] == sim_number + 1) {
      # skip calculation of likelihood for the first observation in each simluation
      sim_number <- sim_number + 1
      next
    }
    
    # observed states and known pars
    x0 <- unlist(data[i-1, "pop"])
    x1 <- unlist(data[i, "pop"])
    sd <- unlist(data[i, sd_est])
    
    # if estimate of sd is invalid skip calculation
    if(is.na(sd) == TRUE || sd <= 0) next
    
    # likelihood calculation
    data[i,"lik0"] <- dnorm(x = x1, mean = x0, sd = sd) # H0: x_t+1 ~ N(x_t, sd)
    data[i,"lik1"] <- likelihood(x_prev = x0, x_now = x1, sd = sd,  # p(x_now | x_prev, disp, theta \in Theta_1)
                                 interval = theta1, guess = guess1)
  }
  
  return(data)
}
