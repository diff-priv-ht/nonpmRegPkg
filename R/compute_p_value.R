#' Compute p-value from regression
#'
#' Function that computes the p-value for a single regression coefficient.
#'
#' @param df A dataframe with the response variable as the first column and the
#'   remaining columns explanatory variables.
#' @param M The number of subsamples
#' @param epsilon The privacy parameter
#' @param beta_number The regression coefficient for which to calculate the
#'   p-value. Numering begins at zero (the intercept).
#' @return The output will be a double
compute_p_value <- function(df, M, epsilon, beta_number){
  n <- nrow(df)
  groups <- sample(rep(1:M, ceiling(n/M))[1:n])
  func <- function(group_num){
    X <- as.matrix(df[groups == group_num,-1])
    Y <- as.matrix(df[groups == group_num, 1])
    betas <- solve(t(X) %*% X) %*% t(X) %*% Y
    betas[beta_number+1]
  }
  coefficients <- map_dbl(.x = 1:M, .f = func)

  p_val <- sum(coefficients > 0) %>%
    Calculate_Z(epsilon = epsilon) %>%
    DP_Binom_test_twoside(n = M, theta_0 = 1/2, epsilon = epsilon, Z = .)
  return(p_val)
}
