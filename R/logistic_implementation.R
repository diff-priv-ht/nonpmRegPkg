#' Test of tests logistic p-value
#'
#' Function that computes the p-value from a single run of the Test of Tests
#' logistic regression.
#'
#' @param df A dataframe with the response variable as the first column and the
#'   remaining columns explanatory variables.
#' @param groups A vector with the group number of each observation.
#' @param M The number of subsamples.
#' @param theta_0 The threshold.
#' @param epsilon The privacy parameter.
#' @return The output will be a double between 0 and 1.
#'
#' @export
logistic_p_val <- function(df, groups, M, theta_0, epsilon){
  func <- function(group_num){
    m <- glm(data = df, Y ~ ., subset = groups == group_num)
    t <- tidy(m)$p.value[2]
    return(t)
  }
  reg_p_vals <- map_dbl(.x = 1:M, .f = func)

  p_val <- sum(reg_p_vals < theta_0, na.rm = TRUE) %>%
    Calculate_Z(epsilon = epsilon) %>%
    DP_Binom_test_oneside(n = M, theta_0 = theta_0, epsilon = epsilon)
  return(2*min(p_val, 1 - p_val))
}


#' Test of Tests Logistic Power
#'
#' Simulate a p-value computation for logistic regression many times and compute
#' the power. The design matrix is regenerated for each simulation.
#'
#' @param nsims Number of simulations.
#' @param n Number of observations (number of rows of the design matrix)
#' @param p Number of parameters (number of columns of the design matrix).
#'   Includes the intercept.
#' @param groups A vector with the group number of each observation.
#' @param epsilon The privacy parameter.
#' @param effect_size The quotient of the parameter of interest (beta) and the
#'   standard deviation of the noise (sigma).
#' @param M The number of subsamples when computing p-value.
#' @param theta_0 The threshold.
#' @param alpha The significance level.
#'
#' @return Will return a vector of p-values.
#'
#' @importFrom parallel mclapply
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @importFrom stats glm
#' @importFrom broom tidy
#'
#' @export
logistic_power <- function(nsims, n, p, groups, epsilon, effect_size, M, theta_0,
                            alpha){
  true_betas <- rep(effect_size, p)

  func <- function(nsim){
    X <- matrix(c(rep(1, n), rnorm(n*(p-1))), ncol = p)
    mean_mat <- X %*% true_betas
    pr = as.vector(1/(1+exp(-mean_mat)))
    Y <- rbinom(nrow(X), size = 1, pr)
    df <- data.frame(Y, X[,-1])
    return(logistic_p_val(df = df, groups = groups, M = M,
                                     epsilon = epsilon, theta_0 = theta_0))
  }
  p_vals <- mclapply(X = 1:nsims, FUN = func, mc.cores = 8)
  return(mean(p_vals < alpha))
}
