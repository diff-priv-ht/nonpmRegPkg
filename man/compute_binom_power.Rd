% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theoretical_power.R
\name{compute_binom_power}
\alias{compute_binom_power}
\title{Power of private binomial test}
\usage{
compute_binom_power(alpha, M, theta_0, thetas, epsilon, nsims, ncores = 1)
}
\arguments{
\item{alpha}{The significance level.}

\item{M}{The number of subsamples.}

\item{theta_0}{The proportion TRUE under the null hypothesis.}

\item{thetas}{A vector of length M with the proportion TRUE in each
subsample.}

\item{epsilon}{The privacy parameter.}

\item{nsims}{The number of draws from the tulap and binomial with which to
compute the reference distribution.}

\item{ncores}{The number of cores to use for the Poisson-binomial pmf
computation}
}
\value{
The output will be a double between 0 and 1.
}
\description{
Function that computes the power of the (two-sided) differentially private
binomial test, adapted from Awan and Slavkovic (2019).
}
