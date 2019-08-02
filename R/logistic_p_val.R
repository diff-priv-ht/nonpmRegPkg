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
#' @importFrom magrittr %>%
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
