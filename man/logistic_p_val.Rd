% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic_implementation.R
\name{logistic_p_val}
\alias{logistic_p_val}
\title{Test of tests logistic p-value}
\usage{
logistic_p_val(df, groups, M, theta_0, epsilon)
}
\arguments{
\item{df}{A dataframe with the response variable as the first column and the
remaining columns explanatory variables.}

\item{groups}{A vector with the group number of each observation.}

\item{M}{The number of subsamples.}

\item{theta_0}{The threshold.}

\item{epsilon}{The privacy parameter.}
}
\value{
The output will be a double between 0 and 1.
}
\description{
Function that computes the p-value from a single run of the Test of Tests
logistic regression.
}
