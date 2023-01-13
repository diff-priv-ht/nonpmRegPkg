#' Function for parameter tuning in Pena-Barrientos
#'
#' The following function is used to find the parameters k, alpha0, and p
#' for the Pena-Barrientos test
#'
#' @param k Determines the number of groups, M, in SSA where M = 2k+1
#' @param p The probability of flipping each bit for the binomial test
#'
#' @importFrom PoissonBinomial ppbinom
#' @importFrom stats pbinom
#'
#' @export
find_eps = function(k, p) {
  t = k
  t_ast = max(t, 2*k-t)
  prob = rep(c(p, 1-p), c(1, 2*k))
  log_pb1 = ppbinom(t_ast, prob, lower.tail = FALSE, log.p = TRUE)
  log_pb0 = pbinom(t_ast, 2*k+1, 1-p, lower.tail = FALSE, log.p = TRUE)
  log_pb1-log_pb0
}

#' Function for parameter tuning in Pena-Barrientos
#'
#' The following function is used to find the parameters k, alpha0, and p
#' for the Pena-Barrientos test
#'
#' @param k Determines the number of groups, M, in SSA where M = 2k+1
#' @param p The probability of flipping each bit for the binomial test
#' @param alpha0 The cuttoff for each test in SSA
#'
#' @export
pow = function(alpha0, k, p) {
  pr = (1-p) + (2*p-1)*alpha0
  pbinom(k, 2*k+1, pr, lower.tail = FALSE)
}

#' Function for parameter tuning in Pena-Barrientos
#'
#' The following function is used to find the parameters k, alpha0, and p
#' for the Pena-Barrientos test
#'
#' @param k Determines the number of groups, M, in SSA where M = 2k+1
#' @param p The probability of flipping each bit for the binomial test
#' @param alpha0 The cuttoff for each test in SSA
#' @param alpha The significance level
#'
#' @export
find_alpha0 = function(alpha0, k, p, alpha) {
  pow(alpha0, k, p)-alpha
}

#' Function for parameter tuning in Pena-Barrientos
#'
#' The following function is used to find the parameters k, alpha0, and p
#' for the Pena-Barrientos test
#'
#' @param k Determines the number of groups, M, in SSA where M = 2k+1
#' @param p The probability of flipping each bit for the binomial test
#' @param eps The privacy parameter
#'
#' @export
f_eps = function(p, k, eps) {
  find_eps(k,p)-eps
}

#' Function for parameter tuning in Pena-Barrientos
#'
#' The following function is used to find the parameters k, alpha0, and p
#' for the Pena-Barrientos test
#'
#' @param k Determines the number of groups, M, in SSA where M = 2k+1
#' @param eps The privacy parameter
#'
#' @importFrom stats uniroot
#'
#' @export
find_p = function(k, eps) {
  tol = 1e-4
  uniroot(f_eps, interval = c(1/2, 1-tol), k = k, eps = eps)$root
}

#' Function for parameter tuning in Pena-Barrientos
#'
#' The following function is used to find the parameters k, alpha0, and p
#' for the Pena-Barrientos test
#'
#' @param alpha The significance level
#' @param eps The privacy parameter
#' @param alpha_min The minimum allowable value for alpha0
#'
#' @export
find_k_alpha0 = function(eps, alpha, alpha_min = 0) {
  k = 0
  p = find_p(k, eps)
  while ( sign(find_alpha0(0, k, p, alpha)) == sign(find_alpha0(1, k, p, alpha)) ) {
    k = k+1
    p = find_p(k, eps)
  }
  alpha0 = uniroot(find_alpha0, interval = c(0,1), alpha = alpha, k = k, p = p)$root
  while (alpha0 < alpha_min) {
    k = k+1
    p = find_p(k, eps)
    alpha0 = uniroot(find_alpha0, interval = c(0,1), alpha = alpha, k = k, p = p)$root
  }
  list(k = k, alpha0 = alpha0)
}

#' Exact power of PB for a normal test
#'
#' The following function gives the exact theoretical power of the Pena-Barrientos
#' method for a normal test.
#'
#' @param epsilon The privacy parameter
#' @param effect_size The quotient of the parameter of interest (mu) and
#'   the standard deviation of the noise (sigma)
#' @param n The sample size (number of rows)
#' @param d The dimension (number of columns)
#' @param n_zeros number of dimensions set to zero
#' @param alpha The significance level
#' @param k,alpha0 Parameters of PB. If NULL, computed within the function
#'
#' @importFrom stats dbinom
#'
#' @export
PB_power_norm <- function(epsilon, effect_size, n, d = 1, n_zeros = 0, alpha = 0.05,
                          k = NA, alpha0 = NA){
  if(is.na(k) | is.na(alpha0)){
    k_alpha0 <- find_k_alpha0(eps = epsilon, alpha = alpha)
    k <- k_alpha0$k; alpha0 <- k_alpha0$alpha0
  }
  p <- find_p(k, epsilon)

  norm_pow <- public_power_normal(n = n/(2*k+1), d = d, effect_size = effect_size,
                                  alpha = alpha0, n_zeros = n_zeros)
  power <- 0
  for(t in 0:(2*k+1)){
    likelihood <- dbinom(t, 2*k+1, norm_pow)
    rej_prob <- ppbinom(x = k, probs = c(rep(p, t), rep(1-p, 2*k+1-t)), lower.tail = F)
    power <- power + likelihood*rej_prob
  }
  return(power)
}

#' Exact power of PB for an ANOVA
#'
#' The following function gives the exact theoretical power of the Pena-Barrientos
#' method for a one-way ANOVA.
#'
#' @param epsilon The privacy parameter
#' @param effect_size The ratio of the between group variance to the within group
#' variance.
#' @param n The sample size (number of rows)
#' @param groups an integer of the number of groups
#' @param alpha The significance level
#' @param k,alpha0 Parameters of PB. If NULL, computed within the function
#'
#' @export
PB_power_ANOVA <- function(epsilon, effect_size, n, groups = 3, alpha = 0.05,
                           k = NA, alpha0 = NA){
  if(is.na(k) | is.na(alpha0)){
    k_alpha0 <- find_k_alpha0(eps = epsilon, alpha = alpha)
    k <- k_alpha0$k; alpha0 <- k_alpha0$alpha0
  }
  p <- find_p(k, epsilon)

  ANOVA_pow <- power.anova.test(groups = groups, between.var = effect_size,
                                within.var = 1, n = floor(n/(2*k+1)/groups),
                                sig.level = alpha0, power = NULL)$power
  power <- 0
  for(t in 0:(2*k+1)){
    likelihood <- dbinom(t, 2*k+1, ANOVA_pow)
    rej_prob <- ppbinom(x = k, probs = c(rep(p, t), rep(1-p, 2*k+1-t)), lower.tail = F)
    power <- power + likelihood*rej_prob
  }
  return(power)
}
