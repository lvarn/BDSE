% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions_processing_results.R
\name{get_counts}
\alias{get_counts}
\title{Extract posterior counts}
\usage{
get_counts(
  results_list,
  rename_cols = FALSE,
  n_burnin = 0,
  n_thin = 1,
  as_array = FALSE
)
}
\arguments{
\item{results_list}{A named list of lists, each elemnt of which corresponds
to a single chain of the Gibbs sampler.}

\item{rename_cols}{A logical scalar indicating whether the count columns
should be renamed (\code{1:n}); defaults to \code{FALSE}. (When
\code{FALSE}, the column names are of the form \code{chainx.counts_xxx},
where \code{x} is a number.)}
}
\value{
A data.frame of posterior count columns.
}
\description{
Extracts the posterior counts vectors for the covariate combinations.
}
