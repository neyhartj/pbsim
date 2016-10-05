#' Define a function to generate random samples from a positive-Poisson distribution
#' 
#' @description This is a hidden utility function that generates samples based on
#' a positive-Poisson distribution. That is, the probability of 0 is 0.
#' 
#' @param n The number of random samples from the distribution
#' @param lambda Mean parameter
#' 
#' @return A numeric vector of random samples from the distribution
#' 
rpospois <- function(n, lambda) {
  qpois(runif(n, min = dpois(0, lambda), max = 1), lambda)
}

