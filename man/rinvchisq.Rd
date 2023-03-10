% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{rinvchisq}
\alias{rinvchisq}
\title{Draw from scaled inverse Chi-squared distribution}
\usage{
rinvchisq(n, deg_free, ssq)
}
\arguments{
\item{n}{An integer specifying the number of observations.}

\item{deg_free}{An integer specifying the degrees of freedom.}

\item{ssq}{The scale parameter, \eqn{s^2}.}
}
\value{
A vector of length \code{n}.
}
\description{
Generates samples from the Scaled Inverse Chi-squared distribution.
}
\details{
Samples are drawn from the Gamma distribution with shape-scale
parameterization:
\deqn{
  Gamma( r / 2, 2 / (r * s^2) ) ,
}
and then transformed. When X has the above Gamma distribution, then
\eqn{1 / X} has a scaled inverse Chi-squared distribution with \eqn{r}
(\code{deg_free}) degrees of freedom and scale parameter \eqn{s^2}
(\code{ssq}).

When \code{ssq} is set to \code{1 / deg_free}, then the drawn samples are
from an inverse Chi-squared distribution with degrees of freedom
\code{deg_free}.
}
