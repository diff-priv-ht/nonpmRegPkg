#' Random draw from Tulap distribution.
#'
#' Created by Canyon.
#'
#' @param n number of observations.
#' @param m,b,q the parameters of the Tulap distribution; described in Awan and
#'   Slavkovic (2018).
#' @return a vector of length \code{n}.
rtulap <- function(n, m, b, q) {
  vect <- rep(NA, n)
  G1 <- rgeom(n, 1 - b)
  G2 <- rgeom(n, 1 - b)
  U <- runif(n, -.5, .5)
  vect <- G1 - G2 + U +m
  return(vect)}

#' CDF of the Tulap distribution.
#'
#' Created by Canyon.
#'
#' @param m,b,q the parameters of the Tulap distribution; described in Awan and
#'   Slavkovic (2018)
tulap_cdf <- function(m, b, x) {
  if (x <= round(m)) {
    (b ^ (- round(x-m))/ (1 + b)) * (b + (x - m - round(x - m) + 1/2) * (1 - b))
  } else {
    1 - ((b ^ round(x - m)) / (1 + b)) * (b + (round(x - m) - (x - m) + 1/2) * (1 - b))  }}

#' Privitize binomial data
#'
#' This function takes in data from a binomial distribution and outputs a
#' differentally private estimate of that value. Created by Canyon.
#'
#' @param X data from a binomial distribution - typically the number of successes.
#' @param epsilon the privacy parameter.
Calculate_Z <- function(X, epsilon) {
  Z <- rtulap(1, X, exp(-epsilon), 0)
  return(Z)}

#' One-sided private binomial hypothesis test.
#'
#' This function takes in a privatized estimate of a draw from a binomial
#' distribution and performs a one-sided hypothesis test. Created by Canyon.
#'
#' @param n number of trials.
#' @param theta_0 null hypothesis.
#' @param epsilon privacy parameter.
#' @param Z privatized number of successes.
#' @return outputs the p-value of the hypothesis test.
DP_Binom_test_oneside <- function(n, theta_0, epsilon, Z) {
  F_underbar <- 0:n %>%
    map_dbl(~ tulap_cdf(0, exp(-epsilon), .x - Z))
  B_underbar <- 0:n %>%
    map_dbl(~ choose(n, .x) * theta_0^.x * (1 - theta_0)^(n - .x))
  return(sum(F_underbar * B_underbar))}

#' Two-sided private binomial hypothesis test.
#'
#' This function takes in a privatized estimate of a draw from a binomial
#' distribution and performs a two-sided hypothesis test. Created by Canyon.
#'
#' @param n number of trials.
#' @param theta_0 null hypothesis.
#' @param epsilon privacy parameter.
#' @param Z privatized number of successes.
#' @return outputs the p-value of the hypothesis test.
DP_Binom_test_twoside <- function(n, theta_0, epsilon, Z) {
  T_stat <- abs(Z - n * theta_0)
  p1 <- DP_Binom_test_oneside(n, theta_0, epsilon, T_stat + n * theta_0)
  p2 <- DP_Binom_test_oneside(n, theta_0, epsilon, - T_stat + n * theta_0)
  return(p1 + 1 - p2)}
