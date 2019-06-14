#' Repeatedly simulate a p-value computation.
#'
#' Simulate a p-value computation many times using choice of computation method
#'
#' @param nsims Number of simulations.
#' @param X Design matrix from which to compute regression coefficients.
#' @param true_betas Vector of true regression parameters. Must be the same
#'   length as the number of columns of \code{X}.
#' @param sigma Standard deviation of the error in the regression.
#' @param M The number of subsamples when computing p-value.
#' @param epsilon The privacy parameter.
#' @param beta_number The regression coefficient for which to calculate the
#'   p-value. Numering begins at zero (the intercept).
#' @param p_value_function Function to use to compute the p-values.
#' @param ... Additional options for \code{p_value_function}.
#'
#' @return Will return a vector of p-values.
#'
#' @export
simulate_p_vals <- function(nsims, X, true_betas, sigma, M, epsilon, beta_number,
                                 p_value_function = compute_p_value, ...){
  n <- nrow(X)
  mean_mat <- X %*% true_betas
  func <- function(){
    Y <- mean_mat + rnorm(n, sd = sigma)
    df <- data.frame(Y, X)
    return(p_value_function(df = df, M = M, epsilon = epsilon, beta_number = beta_number, ...))
  }
  p_vals <- replicate(n = nsims, expr = func())
  return(p_vals)
}
