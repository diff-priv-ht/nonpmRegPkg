#' Power of private binomial test
#'
#' Function that computes the power of the (two-sided) differentially private
#' binomial test, adapted from Awan and Slavkovic (2019).
#'
#' @param alpha The significance level.
#' @param M The number of subsamples.
#' @param theta_0 The proportion TRUE under the null hypothesis.
#' @param thetas A vector of length M with the proportion TRUE in each
#'   subsample.
#' @param epsilon The privacy parameter.
#' @param nsims The number of draws from the tulap and binomial with which to
#'   compute the reference distribution.
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom poibin dpoibin
#' @importFrom stats quantile
#' @importFrom stats rbinom
#'
#' @export
compute_binom_power <- function(alpha, M, theta_0, thetas, epsilon, nsims){
  critical_values <- rbinom(nsims, M, theta_0) %>%
    rtulap(n = nsims, b = exp(-epsilon), q = 0) %>%
    quantile(c(alpha/2, 1- alpha/2))

  cdf_Z <- function(x){
    F_underbar <- 0:M %>%
      map_dbl(~ tulap_cdf(.x, exp(-epsilon), x))
    B_underbar <- 0:M %>%
      map_dbl(~ dpoibin(kk = .x, pp = thetas))
    return(sum(F_underbar * B_underbar))
  }

  return(cdf_Z(as.double(critical_values[1])) +
           (1 - cdf_Z(as.double(critical_values[2]))))
}

#' Theoretical power of our test
#'
#' Function that computes the power of our test for a given a design matrix and
#' a given partitioning into subsamples.
#'
#' @param X For regression only. A design matrix with at least two explanatory
#'   variables.
#' @param groups For regression only. A vector of length \code{nrow(X)} with the
#'   index of the group of each row in \code{X}.
#' @param n For normal test only. The number of observations (number of rows in
#'   the database).
#' @param d For normal test only. The number of dimensions (number of columns in
#'   the database).
#' @param n_zeros For normal test only. The number of entries of the alternative
#'   distribution with mean zero. Defaults to 0.
#' @param M The number of subsamples to partition the data into.
#' @param effect_size The quotient of the parameter of interest (beta) and the
#'   standard deviation of the noise (sigma).
#' @param alpha The significance level.
#' @param epsilon The privacy parameter.
#' @param nsims The number of draws from the tulap and binomial with which to
#'   compute the reference distribution.
#' @param theta_0 The threshold.
#' @param test The test to compute the power of. Either "Linear Regression" or
#'   "Normal"
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#'
#' @export
theoretical_power <- function(X = NULL, groups = NULL, n = NULL, d = NULL,
                              n_zeros = 0, M, effect_size, alpha, epsilon, nsims,
                              theta_0, test = "Linear Regression"){
  if(test == "Linear Regression"){
    X_list <- map(.x = 1:M, .f = function(j){X[groups == j,]})
    thetas <- map_dbl(.x = X_list, .f = public_power_lm, effect_size = effect_size,
                      alpha = theta_0)
  }
  if(test == "Normal"){
    pub_pows <- map_dbl(.x = c(ceiling(n/M), floor(n/M)), .f = public_power_normal,
                        d = d, effect_size = effect_size, alpha = theta_0,
                        n_zeros = n_zeros)
    thetas <- c(rep(pub_pows[1], n - floor(n/M)*M),
                rep(pub_pows[2], M - n + floor(n/M)*M))
  }
  compute_binom_power(alpha = alpha , M = M, thetas = thetas, theta_0 = theta_0,
                      epsilon = epsilon, nsims = nsims)
}
