#' Compute power of private binomial test
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
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#' @importFrom poibin dpoibin
#'
#' @export
compute_binom_power <- function(alpha, M, theta_0, thetas, epsilon, nsims){
  critical_values <- rbinom(nsims, M, theta_0) %>%
    rtulap(nsims, ., exp(-epsilon), 0) %>%
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

#' Compute proportion below the threshold for linear regression
#'
#' Function that computes the proportion of the p-values below a pre-determined
#' threshold, theta_0, for a given subsample of a design matrix in a regression.
#'
#' @param num The index of the subsample of interest.
#' @param X A design matrix for a regression. Must have at least two explanatory
#'   variables.
#' @param groups A vector of length \code{nrow(X)} with the index of the group
#'   of each row in \code{X}.
#' @param effect_size The quotient of the parameter of interest (beta) and the standard
#'   deviation of the noise (sigma).
#' @param theta_0 The threshold
#' @return The output will be a double between 0 and 1.
#'
#' @export
compute_prop_below_threshold_lm <- function(num, X, groups, effect_size, theta_0){
  X_ <- X[groups == num,]
  df <- nrow(X_) - ncol(X_)
  ncp <- effect_size/sqrt(solve(t(X_) %*% X_)[2,2])

  return(pt(qt(theta_0/2, df = df), df = df, ncp = ncp) +
           1 - pt(qt(1 - theta_0/2, df = df), df = df, ncp = ncp))
}

#' Compute proportion below the threshold for the normal test
#'
#' Function that computes the proportion of the p-values below a pre-determined
#' threshold, theta_0, for a given subsample of a design matrix in a test of the
#' mean of multivariate normal data.
#'
#' @param n The number of observations (number of rows in the database).
#' @param d The number of dimensions (number of columns in the database).
#' @param effect_size Determines the mean of the alternate distribution (which
#'   will be `d` repetitions of `effect_size`).
#' @param theta_0 The threshold
#' @return The output will be a double between 0 and 1.
#'
#' @export
compute_prop_below_threshold_normal <- function(n, d, effect_size, theta_0){
  lambda = n*d*effect_size^2

  return(1-pchisq(qchisq(1-theta_0, df = d), df = d, ncp = lambda))
}

#' Compute the theoretical power of our test
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
#' @importFrom magrittr %>%
#'
#' @export
theoretical_power <- function(X = NULL, groups = NULL, n = NULL, d = NULL, M,
                              effect_size, alpha, epsilon, nsims, theta_0,
                              test = "Linear Regression"){
  if(test == "Linear Regression"){
    thetas <- map_dbl(.x = 1:M, .f = compute_prop_below_threshold_lm, X = X,
          groups = groups, effect_size = effect_size, theta_0 = theta_0)
  }
  if(test == "Normal"){
    sizes <- c(rep(ceiling(n/M), n - floor(n/M)*M),
                rep(floor(n/M), M - n + floor(n/M)*M))
    thetas <- map_dbl(.x = sizes, .f = compute_prop_below_threshold_normal,
                      d = d, effect_size = effect_size, theta_0 = theta_0)
  }
  compute_binom_power(alpha = alpha , M = M, thetas = thetas, theta_0 = theta_0,
                      epsilon = epsilon, nsims = nsims)
}
