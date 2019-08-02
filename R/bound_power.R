#'Bound on test's power
#'
#'Function that, for given parameters of the public test, will compute an upper
#'bound on the power of the private test.
#'
#'@param theta,theta_0,n Characteristics of the public test: it requires a
#'  database of size `n` to achieve `theta` power at a sigificance level of
#'  `theta_0`.
#'@param epsilon The privacy parameter.
#'@param alpha The significance level for the private test.
#'@param power The desired power of the private test.
#'@param nsims The number of draws to compute the reference distribution for the
#'  binomial test
#'
#'@export
bound_power <- function(theta, theta_0, epsilon, n, alpha, power, nsims){
  M <- 0
  power_M <- 0
  while(power_M < power){
    M <- M + 1
    power_M <- compute_binom_power(alpha, M, theta_0, thetas = rep(theta, M), epsilon,
                                      nsims = nsims)
  }
  return(paste0("M = ", M, ", n = ", n*M, ", power = ", round(power_M, 2)))
}
