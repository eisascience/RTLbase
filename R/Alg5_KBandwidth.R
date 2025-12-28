

#' Kernel bandwidth estimation for Algorithm 5
#'
#' @description
#' Computes a Silverman-inspired bandwidth to smooth kernel density estimates of
#' projection counts used in Algorithm 5.
#'
#' @param s_k Vector of grid locations.
#' @param c_k Vector of counts or weights at each grid location.
#'
#' @return Numeric bandwidth estimate.
#' @export
KBand_fx <- function(s_k, c_k){
  #this is Algorithm 5, but its faster here as a fx
  #This kernel choice is motivated from the rule of thumb
  #kernel density estimation suggested in Silverman (1986)

  #s_k <- s_j; c_k <- c_j
  if (length(s_k) != length(c_k)) {
    warning("s_k and c_k must have the same length.", call. = FALSE)
    return(NA_real_)
  }
  if (anyNA(s_k) || anyNA(c_k)) {
    warning("s_k and c_k must not contain NA values.", call. = FALSE)
    return(NA_real_)
  }
  if (!length(s_k)) {
    warning("Empty inputs supplied to KBand_fx.", call. = FALSE)
    return(NA_real_)
  }
  N <- sum(c_k)
  s_avg <- sum(s_k*c_k)/N
  sigma_hat <- (sum(c_k*(s_k - s_avg)^2)/(N-1))^0.5
  h <- 0.9*sigma_hat*N^(-0.2)
  bandwidth <- structure(h, class = c("rtl_bandwidth", class(h)),
                         metadata = list(mean = s_avg, sigma = sigma_hat, effective_n = N))
  return(bandwidth)
}
