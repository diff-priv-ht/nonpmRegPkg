#' Barreintos reference distribution
#'
#' Computes the empirical CDF of the reference distribution for the test
#' statistic of the test proposed by Barrientos et al.
#'
#' @param epsilon The privacy parameter.
#' @param n Number of observations in the database (number of rows of the design
#'   matrix).
#' @param a Truncation parameter.
#' @param p Number of parameters (number of columns of the design matrix).
#' @param M Number of subsamples.
#' @param reps Number of repetitions.
#'
#' @return a function with outputs on the interval (0,1)
#'
#' @importFrom dplyr if_else
#' @importFrom rmutil rlaplace
#' @importFrom stats rt
#' @importFrom stats ecdf
#' @importFrom plyr .
#'
#' @export
cdf_Barr <- function(epsilon, n, a, p, M, reps){
  #function to compute the test statistic for a given run
  func <- function(rep){
    coefficients <- rt(M, df = round(n/M) - p) %>%
      if_else(. < a, ., a) %>%
      if_else(. > -a, ., -a)
    sum(coefficients)/sqrt(M) + rlaplace(1, m = 0, s = 2*a/epsilon/sqrt(M))
  }
  #compute the test statistic for reps number of runs and return the empirical cdf
  map_dbl(.x = 1:reps, .f = func) %>%
    abs() %>%
    ecdf()
}

#' Barrientos p-value
#'
#' Function that computes the p-values from a running of the test proposed by
#' Barrientos et al. for several choices of the truncation parameter.
#'
#' @param df A dataframe with the response variable as the first column and the
#'   remaining columns explanatory variables.
#' @param groups A vector with the group number of each observation.
#' @param M The number of subsamples.
#' @param a_vector Vector of values of the truncation parameter.
#' @param epsilon The privacy parameter.
#' @param cdfs List of the CDFs for each level of truncation. Must be of equal
#'   length to `a_vector`.
#' @return The output will be a vector of length `a_vector`, where each entry
#'   correpsonds to the p-value from the test with that level of truncation.
#'
#' @importFrom stats lm
#' @importFrom purrr map2_dbl
#'
#' @export
p_val_Barr <- function(df, groups, M, a_vector, epsilon, cdfs){
  n <- nrow(df)
  #function to run the regression for a given group and return the t-statistic
  func <- function(group_num){
    m <- lm(data = df, Y ~ ., subset = groups == group_num)
    t <- tidy(m)$statistic[2]
  }
  #compute the t-stat for each group, truncate them, and then compute the test stat
  compute_p_val <- function(a, cdf){
    coefficients <- map_dbl(.x = 1:M, .f = func) %>%
      if_else(. < a, ., a) %>%
      if_else(. > -a, ., -a)
    stat <- sum(coefficients)/sqrt(M) + rlaplace(1, m = 0, s = 2*a/epsilon/sqrt(M))

    #return the p-value
    1 - cdf(stat)
  }
  map2_dbl(.x = a_vector, .y = cdfs, .f = compute_p_val)
}

#' Barrientos power
#'
#' Simulate a p-value computation for the test proposed by Barrientos et al.
#' many times and compute the power. This is done for several choices of the
#' truncation parameter, a. The design matrix is set prior to the simulations.
#'
#' @param nsims Number of simulations for the p-value computation.
#' @param reps Number of simulations for the reference distribution computation.
#' @param X The design matrix.
#' @param groups A vector with the group number of each observation.
#' @param epsilon The privacy parameter.
#' @param effect_size The quotient of the parameter of interest (beta) and the
#'   standard deviation of the noise (sigma).
#' @param M The number of subsamples when computing p-value.
#' @param a_vector Vector of truncation levels.
#' @param alpha The significance level.
#'
#' @return Will return a vector of powers of the same length as `a_vector`,
#'   where output corresponds to the power with that level of truncation.
#'
#' @importFrom purrr map
#' @importFrom purrr transpose
#'
#' @export
Barrientos_power <- function(X, groups, M, a_vector, effect_size, alpha, epsilon, nsims,
                          reps){
  n <- nrow(X)
  p <- ncol(X)
  cdfs <- map(.x = a_vector, .f = cdf_Barr, epsilon = epsilon, n = n, p = p, M = M,
              reps = reps)

  mean_mat <- X %*% rep(effect_size, p)
  func <- function(nsim){
    Y <- rnorm(n, mean = mean_mat, sd = 1)
    df <- data.frame(Y, X[,-1])
    return(p_val_Barr(df = df, groups = groups, M = M, a_vector = a_vector,
                      epsilon = epsilon, cdfs = cdfs))
  }
  p_vals <- mclapply(X = 1:nsims, FUN = func) %>%
    transpose()

  power <- function(p_vals){
    mean(p_vals < alpha)
  }
  map_dbl(.x = p_vals, .f = power)
}
