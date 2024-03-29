% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theoretical_power.R
\name{practical_power_norm}
\alias{practical_power_norm}
\title{Find the power of our test for multivariate normal data in a practical setting}
\usage{
practical_power_norm(
  effect_size,
  epsilon,
  n,
  effect_size_grid = seq(0, 1.5, 0.1),
  rho = 0.8,
  alpha = 0.05,
  d = 1,
  n_zeros = 0,
  M_start = min(12, floor(n/3)),
  alpha0_start = 0.2,
  M_max = n
)
}
\arguments{
\item{effect_size}{The quotient of the parameter of interest (mu) and the
standard deviation of the noise (sigma).}

\item{epsilon}{The privacy parameter.}

\item{n}{The number of observations (number of rows in
the database).}

\item{effect_size_grid}{Grid of effect sizes to search over}

\item{rho}{The desired power threshold}

\item{alpha}{The significance level. Defaults to 0.05.}

\item{d}{The number of dimensions (number of columns in
the database). Defaults to 1.}

\item{n_zeros}{The number of entries of the alternative
distribution with mean zero. Defaults to 0.}

\item{M_start}{Starting value of M}

\item{alpha0_start}{Starting value of alpha0}

\item{M_max}{The maximum M to search over.}
}
\description{
Function that finds the parameters, M and theta_0, of our test that yield
high power for a given combination of n and epsilon
in a practical setting where the effect size is not known and then computes
the power at a given effect size.
}
