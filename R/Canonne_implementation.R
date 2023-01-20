#' Canonne Implementation
#'
#' Implementation of Algorithm 3 in Canonne et al. (2019), which tests whether
#' multivariate Gaussian data comes is centered at 0.
#'
#' @param data An `n` by `d` dataframe where `n` is the number of observations
#'   and `d` is the number of dimensions.
#' @param epsilon The privacy parameter.
#' @param delta The privacy parameter.
#' @param alpha As in the paper, this is the effect size the test should be able
#'   to detect. THIS IS NOT THE SIGNIFICANCE LEVEL!
#'
#' @return The result is either to "REJECT" or "FAIL TO REJECT" the null
#'   hypothesis.
#'
#' @importFrom pracma dot
#' @importFrom purrr map_dfc
#' @importFrom purrr map_lgl
#'
#' @export
Canonne_test <- function(data, epsilon, delta, alpha){
  n <- nrow(data)
  d <- ncol(data)

  #Line 1
  if(n < max(25*log(d/delta), 5/epsilon*log(1/delta))){return("REJECT")}

  #Line 2
  m <- function(i){
    c <- sum(data[,i] <= 0)
    return(c/n - .5)
  }

  #Line 3
  z1 <- max(abs(map_dbl(.x = 1:d, .f = m))) + rlaplace(n = 1, s = 1/(epsilon*n))

  #Line 5
  c1 <- sqrt(log(d/delta)/n) + log(1/delta)/(n*epsilon)
  if(z1 > c1){return("REJECT")}

  #Line 6
  B <- 3*sqrt(log(n*d/delta))
  truncate <- function(vect){
    if_else(vect < B, vect, B) %>%
      if_else(. > -B, ., -B)
  }
  data <- map_dfc(.x = data, .f = truncate)

  #Line 9
  X_bars <- map_dbl(.x = data, .f = sum)
  z2 <- max(abs(X_bars)) + rlaplace(n = 1, s = 2*B/epsilon)

  #Line 11
  c2 <- 3*sqrt(2*n*log(n*d/delta)*log(d/delta)) +
    6/epsilon*sqrt(log(n*d/delta))*log(1/delta)
  if(z2 > c2){return("REJECT")}

  #Line 12
  sigma <- B*sqrt(8*d*log(5/(4*delta)))/epsilon
  X_tildes <- X_bars + rnorm(n = d, sd = sigma)

  #Line 13
  Delta_delta <- 144*(d*log(d/delta) + d/(n*epsilon^2)*log(1/delta)^2 +
                        sqrt(n*d*log(d/delta)*log(n/delta)) +
                        sqrt(d)/epsilon*log(1/delta)*
                        sqrt(log(n/delta)))*log(n*d/delta)

  #Line 14
  c3 <- Delta_delta + 36*d/epsilon*log(n*d/delta)*sqrt(log(5/(4*delta))*log(n/delta))

  fun <- function(j){
    pracma::dot(as.vector(data[j,], mode = "numeric"), X_tildes) > c3
  }

  check_condition <- map_lgl(.x = 1:n, .f = fun)

  z3 <- sum(check_condition) %>%
    `+`(rlaplace(n = 1, s = 1/epsilon))

  #Line 15
  if(z3 > log(1/delta)/epsilon){return("REJECT")}

  #Lines 16-20
  X_hat <- data
  problem_rows <- (1:n)[check_condition]
  for(j in problem_rows){
    X_hat[j,] <- rnorm(n = d)
  }

  #Line 21
  T_hat <- sum(map_dbl(.x = X_hat, .f = sum)^2 - n)

  #Line 22
  c4 <- (5*Delta_delta+432*d/epsilon*log(n*d/delta)*
           sqrt(log(n/delta)*log(5/(4*delta))))/epsilon
  z4 <- T_hat + rlaplace(n = 1, s = c4)

  #Line 23
  if(z4 > (n^2*alpha^2)/324){return("REJECT")}

  #Line 24
  return("FAIL TO REJECT")
}

#' Canonne Power
#'
#' Determines the power of the test proposed by Canonne et al. by simulation.
#'
#' @param eff Determines the mean of the alternate distribution (which
#'   will be `d - n_zeros` repetitions of `eff`).
#' @param n Number of observations (number of rows of the database).
#' @param n_zeros Number of dimensions (out of `d`) with no effect.
#' @param d Number of dimensions (number of columns of the database).
#' @param epsilon The privacy parameter.
#' @param delta The privacy parameter.
#' @param alpha As in the paper, this is the effect size the test should be able
#'   to detect. THIS IS NOT THE SIGNIFICANCE LEVEL!
#' @param nsims The number of times to run the test.
#' @param mc.cores Number of cores used for parallelization.
#' @param PC flag for whether or not computation is on a PC
#'
#' @return A double between 0 and 1.
#'
#' @export
Canonne_power <- function(eff, n, n_zeros, d, epsilon, delta, alpha, nsims,
                          mc.cores = detectCores() - 1, PC = FALSE){
  epsilon_scaled <- epsilon/5
  delta_scaled <- delta/17

  func <- function(sim){
    X <- data.frame(matrix(rnorm(n = n*d, mean = c(rep(0, n_zeros),
                                                   rep(eff, d-n_zeros))),
                           nrow = n, byrow = T))
    Canonne_test(data = X, epsilon = epsilon_scaled, delta = delta_scaled,
                 alpha = alpha)
  }
  #X <- map(.x = 1:nsims, .f = func)

  if(!PC){
    results <- mclapply(X = 1:nsims, FUN = func, mc.cores = mc.cores)
  }
  else{
    cl <- makeCluster(mc.cores)
    registerDoParallel(cl)
    clusterExport(cl,list('map_dbl', '%>%', 'rlaplace', 'map_dfc', 'if_else',
                          'map_lgl', 'Canonne_test'))
    results <- parLapply(cl, X = 1:nsims, fun = func)
  }
  return(mean(results == "REJECT"))
}
