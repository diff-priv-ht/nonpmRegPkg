% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PB_implementation.R
\name{pow}
\alias{pow}
\title{Function for parameter tuning in Pena-Barrientos}
\usage{
pow(alpha0, k, p)
}
\arguments{
\item{alpha0}{The cuttoff for each test in SSA}

\item{k}{Determines the number of groups, M, in SSA where M = 2k+1}

\item{p}{The probability of flipping each bit for the binomial test}
}
\description{
The following function is used to find the parameters k, alpha0, and p
for the Pena-Barrientos test
}
