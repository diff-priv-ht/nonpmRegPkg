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
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom poibin dpoibin
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats qexp
#' @importFrom stats optimize
#' @importFrom purrr flatten_dbl
#'
#' @export
compute_binom_power <- function(alpha, M, theta_0, thetas, epsilon){
  cdf_Z <- function(x, thetas){
    F_underbar <- 0:M %>%
      map_dbl(~ tulap_cdf(.x, exp(-epsilon), x))
    B_underbar <- lapply(X = 0:M, FUN = dpoibin, pp = thetas) %>%
      flatten_dbl()
    #map_dbl(~ dpoibin(kk = .x, pp = thetas))
    return(sum(F_underbar * B_underbar))
  }

  tol <- 2*qexp(.99,rate = epsilon)
  critical_value <- optimize(function(x){abs(cdf_Z(x, thetas = rep(theta_0, M)) - 1 + alpha)},
                             interval = c(-tol, M+tol))$minimum


  return(1 - cdf_Z(as.double(critical_value), thetas))
}

#' Theoretical power of our test
#'
#' Function that computes the power of our test for a given a design matrix and
#' a given partitioning into subsamples.
#'
#' @param theta_0 The threshold.
#' @param M The number of subsamples to partition the data into.
#' @param effect_size The quotient of the parameter of interest (beta or mu) and
#'   the standard deviation of the noise (sigma). For ANOVA, a vector with the
#'   between group variance and within group variance.
#' @param epsilon The privacy parameter.
#' @param alpha The significance level, defaults to 0.05
#' @param X For regression only. A design matrix with at least two explanatory
#'   variables.
#' @param groups For regression, a vector of length \code{nrow(X)} with the
#'   index of the group of each row in \code{X}. For ANOVA, an integer of the
#'   number of groups
#' @param n For normal or ANOVA. The number of observations (number of rows in
#'   the database).
#' @param d For normal test only. The number of dimensions (number of columns in
#'   the database).
#' @param n_zeros For normal test only. The number of entries of the alternative
#'   distribution with mean zero. Defaults to 0.
#' @param nsims The number of draws from the tulap and binomial with which to
#'   compute the reference distribution. (No Longer Used)
#' @param test The test to compute the power of. Either "Linear Regression",
#'   "Normal", or "ANOVA"
#' @param ncores The number of cores to use for the Poisson-binomial pmf
#'   computation (No Longer Used)
#' @return The output will be a double between 0 and 1.
#'
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom stats power.anova.test
#'
#' @export
theoretical_power <- function(theta_0, M, effect_size, epsilon, alpha = 0.05,
                              X = NULL, groups = NULL, n = NULL, d = 1,
                              n_zeros = 0, nsims = NULL,
                              test = "Linear Regression", ncores = 1){
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
  if(test == "ANOVA"){
    g <- groups
    if(length(effect_size) == 1){ effect_size <- rep(effect_size, 2) }
    f <- function(n){power.anova.test(groups = g, between.var = effect_size[1],
                                      within.var = effect_size[2], n = n,
                                      sig.level = theta_0, power = NULL)$power}
    pub_pows <- map_dbl(.x = c(ceiling(n/M/g), floor(n/M/g)), .f = f)
    thetas <- c(rep(pub_pows[1], n - floor(n/M)*M),
                rep(pub_pows[2], M - n + floor(n/M)*M))
  }
  compute_binom_power(alpha = alpha , M = M, thetas = thetas, theta_0 = theta_0,
                      epsilon = epsilon)
}

#' Optimize the parameters of our test for multivariate normal data
#'
#' Function that finds the parameters, M and theta_0, of our test that yield
#' the highest power for a given combination of n, the effect size, and epsilon
#'
#' @param effect_size The quotient of the parameter of interest (mu) and the
#'   standard deviation of the noise (sigma).
#' @param epsilon The privacy parameter.
#' @param n The number of observations (number of rows in
#'   the database).
#' @param alpha The significance level. Defaults to 0.05.
#' @param d The number of dimensions (number of columns in
#'   the database). Defaults to 1.
#' @param n_zeros The number of entries of the alternative
#'   distribution with mean zero. Defaults to 0.
#' @param M_max The maximum M to search over. The optimizer will try all
#'   integers from 1 to M_max
#' @return The output will be a list with the achievable power and the
#'   corresponding two parameter values
#'
#' @importFrom stats optimize
#'
#' @export
optimize_power_norm <- function(effect_size, epsilon, n, alpha = 0.05, d = 1,
                                n_zeros = 0, M_max = min(floor(n/3), 50)){
  curr_max_pow <- 0

  for(i in 1:M_max){
    opt <- optimize(f = theoretical_power, interval = c(0,1), maximum = T, M = i,
                    effect_size = effect_size, alpha = alpha, epsilon = epsilon,
                    test = "Normal", n = n, d = d, X = NULL, groups = NULL,
                    n_zeros = n_zeros, nsims = NULL)
    if(opt$objective > curr_max_pow){
      curr_max_pow <- opt$objective
      theta_0 <- opt$maximum
      M <- i
    }
  }
  return(list("power" = curr_max_pow, "theta_0" = theta_0, "M" = M))
}


#' Optimize the parameters of our test for multivariate normal data
#'
#' Function that finds the parameters, M and theta_0, of our test that yield
#' the highest power for a given combination of n, the effect size, and epsilon
#'
#' @param n The number of observations (number of rows in the database). Should
#'   be a multiple of three.
#' @param effect_size A vector with the  between group variance and within group
#'   variance. Defaults to 1, as in Couch et al.
#' @param epsilon The privacy parameter. Defaults to 1, as in Couch et al.
#' @param alpha The significance level. Defaults to 0.05.
#' @param groups The number of groups in each subsample for ANOVA.
#' @param M_max The maximum M to search over. The optimizer will try all
#'   integers from 1 to M_max
#' @return The output will be a list with the achievable power and the
#'   corresponding two parameter values
#'
#' @importFrom stats optimize
#'
#' @export
optimize_power_ANOVA <- function(n, effect_size = 1, epsilon = 1, alpha = 0.05,
                                 groups = 3, M_max = min(floor(n/groups/2), 50)){
  curr_max_pow <- 0

  for(i in 1:M_max){
    opt <- optimize(f = theoretical_power, interval = c(0,1), maximum = T, M = i,
                    effect_size = effect_size, alpha = alpha, epsilon = epsilon,
                    test = "ANOVA", n = n, groups = groups, X = NULL,
                    n_zeros = NULL, nsims = NULL)
    if(opt$objective > curr_max_pow){
      curr_max_pow <- opt$objective
      theta_0 <- opt$maximum
      M <- i
    }
  }
  return(list("power" = curr_max_pow, "theta_0" = theta_0, "M" = M))
}
