% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_synthetic_dataset.R
\name{make_observed_dataset}
\alias{make_observed_dataset}
\title{Make observed dataset}
\usage{
make_observed_dataset(target_dataset, clustering_var, clusters_col = "cluster")
}
\arguments{
\item{target_dataset}{A data.frame, the target dataset with cell indicator
and cluster columns.}

\item{clustering_var}{A string specifying the clustering variable (the name
of a column in target dataset).}

\item{clusters_col}{A string specifying the sequential cluster IDs column in
the target dataset. (Defaults to "cluster".)}
}
\value{
A data.frame, the observed dataset. If some clusters are not
  observed, then an additional sequential cluster ID column called
  \code{paste0(clusters_col, "_obs")} will be added.
}
\description{
Makes the observed dataset for the Gibbs sampler by excluding records not on
either of the two lists and creating a sequential clusters column if some
clusters in the target dataset are not observed.
}
