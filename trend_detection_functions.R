# linear trend detection

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

### simulation functions

# linear population trend model w/ normal process error
sim_lin_proc <- function (x0, t, trend, sd, runstat = TRUE, proc_err = TRUE) {
  # allow for constant trend or vector of varying trends (automatically recycles vector if lenght(trend) < t)
  trend_val = rep(trend, length = t)
  sd_val = rep(sd, length = t)
  
  # run model with proc error or w/o proc error
  if (proc_err) {
    pop <- vector("numeric", length = t)
    pop[1] <- x0
    for (i in 1:(length(pop) - 1)) {
      pop[i + 1] <- rnorm(1, pop[i] + trend_val[i], sd_val[i])
    }
  } else {
    pop <- x0 + cumsum(trend_val) + rnorm(t, mean = 0, sd = sd_val)
  }
  
  # create df of results
  out <- tibble(time = 1:t, pop = pop, trend = trend_val, sd = sd_val)
  
  # compute running statistics to use for sequential test calculations
  if(runstat == TRUE){
    out$diff_1 <- out$pop - lag(out$pop)
    running_stats <- run_stat(diff(out$pop))
    out <- dplyr::left_join(out, running_stats, by = "diff_1")
  }
  return(out)
  #return(tibble(time = 1:t, pop = pop, trend = trend_val, sd = sd_val))
}

# add normal and unbiased observation noise to a population time series
obs_model <- function (x, sd) {
  obs <- x + rnorm(n = length(x), mean = 0, sd = sd)
  return(obs)
}

# simulate sim_lin_proc nsim times using common parameters
multi_sim <- function (nsim, x0, t, trend, sd, runstat = TRUE, proc_err = TRUE) {
  simulations <- lapply(1:nsim, function(x) sim_lin_proc(x0, t, trend, sd, runstat, proc_err))
  return(bind_rows(simulations, .id = "simulation") %>% mutate(simulation = as.numeric(simulation)))
}

###
sim_log_gomp <- function (xinit,     # initial log population size
                     lambda, # population growth rate,
                     trend = 0,  # log-linear trend term
                     b,      # density dependence
                     phi,    # process noise autocorrelation [-1,1]
                     sd_proc, # process noise: N(0, sd_proc)
                     sd_obs,  # observation noise: N(0, sd_obs)
                     tfinal,
                     nsim,
                     runstat = TRUE) {
  
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
      log_x[j] <- lambda + b * log_x[j-1] + trend + proc_error
      log_y[j] <- log_x[j] + rnorm(n = 1, mean = 0, sd = sd_obs)
    }
    # save current simulation
    sim_list[[i]] <- (tidyr::tibble(time = 1:tfinal,
                                    log_x = log_x,
                                    log_y = log_y,
                                    lambda = lambda,
                                    trend = trend,
                                    b_ = b,
                                    phi = phi,
                                    sd_proc = sd_proc,
                                    sd_obs = sd_obs))
  } 
  out <- dplyr::bind_rows(sim_list, .id = "simulation")
  
  if(runstat == TRUE){
    out$diff_1 <- out$log_x - lag(out$log_x)
    running_stats <- run_stat(diff(out$log_x))
    out <- dplyr::left_join(out, running_stats, by = "diff_1")
  }
  
  return(out)
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
    data$pval <- NA_real_
    data$trend_est <- NA_real_
    
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

get_conclusion <- function(data, sig_level) {
  data %>% 
    mutate(significant = ifelse(pval <= sig_level, TRUE, FALSE),
           correct_sign = ifelse(sign(trend_est) == sign(trend), TRUE, FALSE),
           correct_inference = ifelse(significant == TRUE & correct_sign == TRUE, TRUE, FALSE),
           type_error = case_when(correct_sign == FALSE & significant == FALSE ~ "insignificant_false",
                                  correct_sign == FALSE & significant == TRUE ~ "significant_false",
                                  correct_sign == TRUE & significant == FALSE ~ "insignificant_true",
                                  correct_sign == TRUE & significant == TRUE ~ "significant_true")) %>%
    return()
}

get_power <- function(data, power_level = 0.8, nsim) {
  data %>%
    group_by(time) %>%
    count(type_error) %>%
    filter(type_error == "significant_true") %>%
    ungroup() %>%
    mutate(time_power = ifelse(n >= nsim * power_level, time, NA),
           mintime = min(time_power, na.rm = TRUE)) %>%
    select(-time_power) %>%
    return()
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
calc_sr <- function(data, accept_reg = c(0,0), reject_reg = c(0, 5), sd_est, sd_min = 0) {
  
  # initialize empty columns for likelihood calculations
  data$lik0 <- NA_real_
  data$lik1 <- NA_real_
  
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
    # replace sd if lower than minimum (default sd_min == 0)
    sd <- max(sd_min, sd)
    
    # likelihood calculation
    data[i,"lik0"] <- dnorm(x = x1, mean = x0, sd = sd) # H0: x_t+1 ~ N(x_t, sd)
    data[i,"lik1"] <- likelihood(x_prev = x0, x_now = x1, sd = sd,  # p(x_now | x_prev, disp, theta \in Theta_1)
                                 interval = theta1, guess = guess1)
  }
  
  return(data)
}

# calc_sprt <- function (data, accept_reg = c(0,0), reject_reg = c(0, 5), sd_est, alpha = 0.05, beta = 0.2) {
#   
#   OUT <- data %>%
#     calc_sr(sd_est = sd_est) %>%
#     group_by(simulation) %>%
#     mutate(sprt = log(lik1/lik0) ,
#            sprt = ifelse(is.na(sprt), 0, sprt),
#            sprt = cumSkipNA(sprt, cum),
#            a = log( (1 - beta)/ alpha),
#            b = log(beta/(1 - alpha)),
#            decision = case_when(sprt <= a & sprt >= b ~"?",
#                                 sprt > a ~ "H1",
#                                 sprt < b ~ "H0"))
#   
#   # # calculate likelihoods under each hypothesis
#   # OUT <- calc_sr(data, accept_reg = c(0,0), reject_reg = c(0, 5), sd_est)
#   # # calculate sequential test statistic
#   # OUT$sprt <- log(OUT$lik1/OUT$lik0)
#   # OUT$sprt <- ifelse(is.na(OUT$sr), 0, OUT$sr) # replace NA with 0 so that I can use cumsum() to get running LR
#   # OUT$sprt <- cumsum(OUT$sprt)
#   # 
#   # # calculate decision thresholds
#   # OUT$a <- rep(log( (1 - beta)/ alpha), length = nrow(OUT))
#   # OUT$b <- rep(log(beta/(1 - alpha)), length = nrow(OUT))
#   # 
#   # # assign decision according to whether SPRT crosses a threshold
#   # OUT$decision <- ifelse(OUT$sprt <= a & OUT$sprt >= b,
#   #                        "?",
#   #                        ifelse(OUT$sprt > a, "H1", "H0"))
#   return(OUT)
# }

calc_sprt <- function (data, accept_reg = c(0,0), reject_reg = c(0, 5), sd_est, alpha = 0.05, beta = 0.2, sd_min = 0) {
  # this function calculates the log-likelhiood ratio as log(p1) - log(p0) instead of using division
  # 2nd step of SPRT calculation now sets all NA values to 0 (we can't compare models when we don't have sd_mle yet)
  # 3rd step of SPRT calculation sets all +/-Inf values to NA. When taking cumulative LR we skip NA values
  OUT <- data %>%
    calc_sr(sd_est = sd_est, sd_min = sd_min) %>%
    group_by(simulation) %>%
    mutate(sprt = log(lik1) - log(lik0) ,
           sprt = ifelse(is.na(sprt), 0, sprt),
           sprt = ifelse(sprt %in% c(-Inf,Inf), NA, sprt),
           sprt = cumSkipNA(sprt, sum),
           a = log( (1 - beta)/ alpha),
           b = log(beta/(1 - alpha)),
           decision = case_when(sprt <= a & sprt >= b ~"?",
                                sprt > a ~ "H1",
                                sprt < b ~ "H0"),
           decision = ifelse(is.na(decision), "?", decision))
  
  # # calculate likelihoods under each hypothesis
  # OUT <- calc_sr(data, accept_reg = c(0,0), reject_reg = c(0, 5), sd_est)
  # # calculate sequential test statistic
  # OUT$sprt <- log(OUT$lik1/OUT$lik0)
  # OUT$sprt <- ifelse(is.na(OUT$sr), 0, OUT$sr) # replace NA with 0 so that I can use cumsum() to get running LR
  # OUT$sprt <- cumsum(OUT$sprt)
  # 
  # # calculate decision thresholds
  # OUT$a <- rep(log( (1 - beta)/ alpha), length = nrow(OUT))
  # OUT$b <- rep(log(beta/(1 - alpha)), length = nrow(OUT))
  # 
  # # assign decision according to whether SPRT crosses a threshold
  # OUT$decision <- ifelse(OUT$sprt <= a & OUT$sprt >= b,
  #                        "?",
  #                        ifelse(OUT$sprt > a, "H1", "H0"))
  return(OUT)
}

# helper function for calculating time to decision
is_h1 <- function(x) ifelse(x == "H1", TRUE, FALSE)
is_h0 <- function(x) x == "H0"

detection_time_est <- function (data, FUNC = is_h1) {
  data %>%
    group_by(simulation) %>%
    summarize(delay_est = unique(detect_index(decision_est, FUNC)),
              delay_est = ifelse(delay_est == 0, NA, delay_est)) %>%
    # mutate(mean_delay_est = mean(delay_est),
    #        median_delay_est = median(delay_est)) %>% #,
    #        p20_delay = quantile(delay, probs = 0.2),
    #        p80_delay = quantile(delay, probs = 0.8),
    #        p90_delay = quantile(delay, probs = 0.9)) %>%
    return()
}

detection_time_known <- function (data, FUNC = is_h1) {
  data %>%
    group_by(simulation) %>%
    summarize(delay_known = unique(detect_index(decision_known, FUNC)),
              delay_known = ifelse(delay_known == 0, NA, delay_known)) %>%
    # mutate(mean_delay_known = mean(delay_known),
    #        median_delay_known = median(delay_known)) %>% #,
    #        p20_delay = quantile(delay, probs = 0.2),
    #        p80_delay = quantile(delay, probs = 0.8),
    #        p90_delay = quantile(delay, probs = 0.9)) %>%
    return()
}
