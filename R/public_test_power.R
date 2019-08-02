#' Linear regression public power
#'
#' Function that computes the power of public linear regression for a given
#' design matrix and effect size.
#'
#' @param X A design matrix for a regression. Must have at least two columns.
#' @param effect_size The quotient of the parameter of interest (beta) and the
#'   standard deviation of the noise (sigma).
#' @param alpha The significance level.
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom stats qt
#' @importFrom stats pt
#'
#' @export
public_power_lm <- function(X, effect_size, alpha){
  df <- nrow(X) - ncol(X)
  ncp <- effect_size/sqrt(solve(t(X) %*% X)[2,2])

  return(pt(qt(alpha/2, df = df), df = df, ncp = ncp) +
           1 - pt(qt(1 - alpha/2, df = df), df = df, ncp = ncp))
}

#' Normal test public power
#'
#' Function that computes the publuc power for the test of the mean of
#' multivariate normal data.
#'
#' @param n The number of observations (number of rows in the database).
#' @param d The number of dimensions (number of columns in the database).
#' @param effect_size Determines the mean of the alternate distribution (which
#'   will be `d` repetitions of `effect_size`).
#' @param alpha The significance level.
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom stats qchisq
#' @importFrom stats pchisq
#'
#' @export
public_power_normal <- function(n, d, effect_size, alpha){
  lambda = n*d*effect_size^2

  return(1-pchisq(qchisq(1-alpha, df = d), df = d, ncp = lambda))
}
