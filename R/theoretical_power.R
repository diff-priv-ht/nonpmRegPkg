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
#'   the standard deviation of the noise (sigma). For ANOVA, the ratio of the
#'   between group variance to the within group variance.
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
    f <- function(n){power.anova.test(groups = g, between.var = effect_size,
                                      within.var = 1, n = n,
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
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel detectCores
#'
#' @export
optimize_power_norm <- function(effect_size, epsilon, n, alpha = 0.05, d = 1,
                                n_zeros = 0, M_max = min(floor(n/3), 50)){

  f <- function(M, n, d, n_zeros, effect_size, epsilon, alpha){
    opt <- optimize(f = theoretical_power, interval = c(0,1), maximum = T, M = M,
                    effect_size = effect_size, alpha = alpha, epsilon = epsilon,
                    test = "Normal", n = n, d = d, X = NULL, groups = NULL,
                    n_zeros = n_zeros, nsims = NULL, tol = 5e-3)
    return(c(opt$objective, opt$maximum,M))
  }

  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  clusterExport(cl,list('theoretical_power', 'map_dbl', 'compute_binom_power',
                        'public_power_normal', '%>%', 'tulap_cdf',
                        'flatten_dbl', 'dpoibin'))
  x <- parLapply(cl, X = 1:M_max, fun = f, n = n, d = d, n_zeros = n_zeros,
                 effect_size = effect_size, epsilon = epsilon, alpha = alpha)
  max_x <- x[[which.max(as.numeric(sapply(x,"[[",1)))]]

  return(list("power" = max_x[1], "theta_0" = max_x[2], "M" = as.integer(max_x[3])))
}

#' Find practical values of the parameters for a normal test
#'
#' Performs a binary search to find the parameter values that reach power rho at
#' the smallest effect size. If no effect sizes are sufficient to reach power
#' rho, returns the parameters for the largest effect size in input grid.
#'
#' @param effect_size_grid Grid of effect sizes to search over
#' @param epsilon The privacy parameter.
#' @param n The number of observations (number of rows in
#'   the database).
#' @param rho The desired power threshold
#' @param alpha The significance level. Defaults to 0.05.
#' @param d The number of dimensions (number of columns in
#'   the database). Defaults to 1.
#' @param n_zeros The number of entries of the alternative
#'   distribution with mean zero. Defaults to 0.
#' @param M_start Starting value of M
#' @param alpha0_start Starting value of alpha0
#'
#' @return The output will be a list with the achievable power, the
#'   corresponding two parameter values, and the minimum effect size
#'
#' @export
practical_pars_norm <- function(effect_size_grid, epsilon, n, rho = 0.8, alpha = 0.05,
                                  d = 1, n_zeros = 0, M_start = min(12,floor(n/3)),
                                  alpha0_start = 0.2){
  i <- (length(effect_size_grid)+1) %/% 2; M = M_start; alpha_0 = alpha0_start
  while(length(effect_size_grid) > 1){
    pwr <- theoretical_power(effect_size = effect_size_grid[i], M = M, n = n,
                             epsilon = epsilon, test = "Normal", theta_0 = alpha_0,
                             d = d, n_zeros = n_zeros, alpha = alpha)
    if(pwr < rho){
      opt <- optimize_power_norm(effect_size = effect_size_grid[i], epsilon, n,
                                 alpha, d, n_zeros)
      pwr <- opt$power; alpha_0 <- opt$theta_0; M <- opt$M
    }
    if(pwr >= rho){
      effect_size_grid <- effect_size_grid[1:i]
      i <- (i+1) %/% 2
    }
    else{
      effect_size_grid <- effect_size_grid[(i+1):length(effect_size_grid)]
      i <- (length(effect_size_grid) + 1) %/% 2
    }
  }
  opt <- optimize_power_norm(effect_size = effect_size_grid, epsilon, n,
                             alpha, d, n_zeros)
  return(list("power" = opt$power, "theta_0" = opt$theta_0, "M" = opt$M,
              "min_effect_size" = effect_size_grid))
}

#' Find the power of our test for multivariate normal data in a practical setting
#'
#' Function that finds the parameters, M and theta_0, of our test that yield
#' high power for a given combination of n and epsilon
#' in a practical setting where the effect size is not known and then computes
#' the power at a given effect size.
#'
#' @param effect_size The quotient of the parameter of interest (mu) and the
#'   standard deviation of the noise (sigma).
#' @param epsilon The privacy parameter.
#' @param n The number of observations (number of rows in
#'   the database).
#' @param effect_size_grid Grid of effect sizes to search over
#' @param rho The desired power threshold
#' @param alpha The significance level. Defaults to 0.05.
#' @param d The number of dimensions (number of columns in
#'   the database). Defaults to 1.
#' @param n_zeros The number of entries of the alternative
#'   distribution with mean zero. Defaults to 0.
#' @param M_start Starting value of M
#' @param alpha0_start Starting value of alpha0
#'
#' @export
practical_power_norm <- function(effect_size, epsilon, n, effect_size_grid = seq(0,1.5,0.1),
                                 rho = 0.8, alpha = 0.05, d = 1, n_zeros = 0,
                                 M_start = min(12,floor(n/3)), alpha0_start = 0.2){
  l <- practical_pars_norm(effect_size_grid, epsilon, n, rho, alpha, d, n_zeros,
                            M_start, alpha0_start)
  return(theoretical_power(effect_size = effect_size, M = l$M, n = n,
                           epsilon = epsilon, test = "Normal", theta_0 = l$theta_0,
                           d = d, n_zeros = n_zeros, alpha = alpha))
}


#' Optimize the parameters of our test for ANOVA data
#'
#' Function that finds the parameters, M and theta_0, of our test that yield
#' the highest power for a given combination of n, the effect size, and epsilon
#'
#' @param n The number of observations (number of rows in the database). Should
#'   be a multiple of three.
#' @param effect_size The ratio of the between group variance to the within group
#' variance. Defaults to 1, as in Couch et al.
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
  f <- function(M, n, groups, effect_size, epsilon, alpha){
    opt <- optimize(f = theoretical_power, interval = c(0,1), maximum = T, M = M,
                    effect_size = effect_size, alpha = alpha, epsilon = epsilon,
                    test = "ANOVA", n = n, groups = groups, X = NULL,
                    n_zeros = NULL, nsims = NULL, tol = 5e-3)
    return(c(opt$objective, opt$maximum,M))
  }

  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  clusterExport(cl,list('theoretical_power', 'map_dbl', 'compute_binom_power',
                        '%>%', 'tulap_cdf', 'flatten_dbl', 'dpoibin'))
  x <- parLapply(cl, X = 1:M_max, fun = f, n = n, groups = groups,
                 effect_size = effect_size, epsilon = epsilon, alpha = alpha)
  max_x <- x[[which.max(as.numeric(sapply(x,"[[",1)))]]

  return(list("power" = max_x[1], "theta_0" = max_x[2], "M" = as.integer(max_x[3])))
}

#' Find practical values of the parameters for an ANOVA
#'
#' Performs a binary search to find the parameter values that reach power rho at
#' the smallest effect size. If no effect sizes are sufficient to reach power
#' rho, returns the parameters for the largest effect size in input grid.
#'
#' @param effect_size_grid Grid of effect sizes to search over
#' @param epsilon The privacy parameter.
#' @param n The number of observations (number of rows in
#'   the database).
#' @param rho The desired power threshold
#' @param alpha The significance level. Defaults to 0.05.
#' @param groups The number of groups in each subsample for ANOVA.
#' @param M_start Starting value of M
#' @param alpha0_start Starting value of alpha0
#'
#' @return The output will be a list with the achievable power, the
#'   corresponding two parameter values, and the minimum effect size
#'
#' @export
practical_pars_ANOVA <- function(effect_size_grid, epsilon, n, rho = 0.8,
                                 alpha = 0.05, groups = 3,
                                 M_start = min(12,floor(n/2/groups)),
                                 alpha0_start = 0.2){
  i <- (length(effect_size_grid)+1) %/% 2; M = M_start; alpha_0 = alpha0_start
  while(length(effect_size_grid) > 1){
    pwr <- theoretical_power(effect_size = effect_size_grid[i], M = M, n = n,
                             epsilon = epsilon, test = "ANOVA", theta_0 = alpha_0,
                             groups = groups, alpha = alpha)
    if(pwr < rho){
      opt <- optimize_power_ANOVA(n, effect_size = effect_size_grid[i], epsilon,
                                 alpha, groups)
      pwr <- opt$power; alpha_0 <- opt$theta_0; M <- opt$M
    }
    if(pwr >= rho){
      effect_size_grid <- effect_size_grid[1:i]
      i <- (i+1) %/% 2
    }
    else{
      effect_size_grid <- effect_size_grid[(i+1):length(effect_size_grid)]
      i <- (length(effect_size_grid) + 1) %/% 2
    }
  }
  if(length(effect_size_grid) == 0){
    print("No Solutions on this grid")
    return(list("power" = NA, "theta_0" = NA, "M" = NA, "min_effect_size" = NA))
  }
  opt <- optimize_power_ANOVA(n, effect_size = effect_size_grid, epsilon,
                              alpha, groups)
  return(list("power" = opt$power, "theta_0" = opt$theta_0, "M" = opt$M,
              "min_effect_size" = effect_size_grid))
}


#' Find the power of our test for an ANOVA in a practical setting
#'
#' Function that finds the parameters, M and theta_0, of our test that yield
#' high power for a given combination of n and epsilon
#' in a practical setting where the effect size is not known and then computes
#' the power at a given effect size.
#'
#' @param effect_size The ratio of the between group variance to the within group
#' variance.
#' @param epsilon The privacy parameter.
#' @param n The number of observations (number of rows in
#'   the database).
#' @param effect_size_grid Grid of effect sizes to search over
#' @param rho The desired power threshold
#' @param alpha The significance level. Defaults to 0.05.
#' @param groups The number of groups in each subsample for ANOVA.
#' @param M_start Starting value of M
#' @param alpha0_start Starting value of alpha0
#'
#' @export
practical_power_ANOVA <- function(effect_size, epsilon, n, effect_size_grid = seq(0,3,0.2),
                                  rho = 0.8, alpha = 0.05, groups = 3,
                                  M_start = min(12,floor(n/2/groups)), alpha0_start = 0.2){
  l <- practical_pars_ANOVA(effect_size_grid, epsilon, n, rho, alpha, groups,
                            M_start, alpha0_start)
  return(theoretical_power(effect_size = effect_size, M = l$M, n = n,
                           epsilon = epsilon, test = "ANOVA", theta_0 = l$theta_0,
                           groups = groups, alpha = alpha))
}
