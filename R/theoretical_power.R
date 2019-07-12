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

#' Compute proportion below the threshold
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
compute_prop_below_threshold <- function(num, X, groups, effect_size, theta_0){
  X_ <- X[groups == num,]
  df <- nrow(X_) - ncol(X_)
  ncp <- effect_size/sqrt(solve(t(X_) %*% X_)[2,2])

  return(pt(qt(theta_0/2, df = df), df = df, ncp = ncp) +
           1 - pt(qt(1 - theta_0/2, df = df), df = df, ncp = ncp))
}

#' Compute the theoretical power of our test
#'
#' Function that computes the power of our test for a given a design matrix and
#' a given partitioning into subsamples.
#'
#' @param X A design matrix for a regression. Must have at least two explanatory
#'   variables.
#' @param groups A vector of length \code{nrow(X)} with the index of the group
#'   of each row in \code{X}.
#' @param M The number of subsamples to partition the data into.
#' @param effect_size The quotient of the parameter of interest (beta) and the
#'   standard deviation of the noise (sigma).
#' @param alpha The significance level.
#' @param epsilon The privacy parameter.
#' @param nsims The number of draws from the tulap and binomial with which to
#'   compute the reference distribution.
#' @param theta_0 The threshold.
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#'
#' @export
theoretical_power <- function(X, groups, M, effect_size, alpha, epsilon, nsims,
                                theta_0){
  map_dbl(.x = 1:M, .f = compute_prop_below_threshold, X = X,
          groups = groups, effect_size = effect_size, theta_0 = theta_0) %>%
    compute_binom_power(alpha = alpha , M = M, thetas = .,
                           theta_0 = theta_0, epsilon = epsilon, nsims = nsims)
}
