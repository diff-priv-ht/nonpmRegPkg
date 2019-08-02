# bound_power <- function(theta, theta_0, epsilon, n, alpha, power, nsims){
#   M <- 0
#   power_M <- 0
#   while(power_M < power){
#     M <- M + 1
#     power_M <- compute_binom_power(alpha, M, theta_0, thetas = rep(theta, M), epsilon,
#                                       nsims = nsims)
#   }
#   return(paste0("M = ", M, ", n = ", n*M, ", power = ", round(power_M, 2)))
# }
#
# bound_power(nsims = 10000, theta = .95, theta_0 = 0.05, epsilon = .1, n = 150, alpha = 0.05, power = .95)
