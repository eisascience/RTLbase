

KBand_fx <- function(s_k, c_k){
  #this is Algorithm 5, but its faster here as a fx
  #This kernel choice is motivated from the rule of thumb
  #kernel density estimation suggested in Silverman (1986)

  #s_k <- s_j; c_k <- c_j
  N <- sum(c_k)
  s_avg <- sum(s_k*c_k)/N
  sigma_hat <- (sum(c_k*(s_k - s_avg)^2)/(N-1))^0.5
  h <- 0.9*sigma_hat*N^(-0.2)
  return(h)
}
