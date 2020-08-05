# linear trend detection

### helper functions

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
  return(tidyr::tibble(x, mu_mle = m_, sd_mle, sd_corrected = sd_) %>% tibble::rowid_to_column(var = "id"))
}

#out %>% ggplot() + geom_line(aes(x = id, y = mu_mle))
run_stat(rnorm(1000))
