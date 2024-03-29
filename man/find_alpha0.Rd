% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PB_implementation.R
\name{find_alpha0}
\alias{find_alpha0}
\title{Function for parameter tuning in Pena-Barrientos}
\usage{
find_alpha0(alpha0, k, p, alpha)
}
\arguments{
\item{alpha0}{The cuttoff for each test in SSA}

\item{k}{Determines the number of groups, M, in SSA where M = 2k+1}

\item{p}{The probability of flipping each bit for the binomial test}

\item{alpha}{The significance level}
}
\description{
The following function is used to find the parameters k, alpha0, and p
for the Pena-Barrientos test
}
