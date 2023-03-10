% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions_synthetic_data.R
\name{get_dummified_names}
\alias{get_dummified_names}
\title{Get names of dummified data.frame}
\usage{
get_dummified_names(dataset, model_formula = ~.)
}
\arguments{
\item{dataset}{A data.frame.}

\item{model_formula}{A formula, which will be passed to \code{model.matrix()}
in the function to dummify the input data.frame. (Defaults to \code{~ .} if
not specified.)}
}
\value{
A single-column character matrix of dummified variable names.
}
\description{
This function returns the names of the input data.frame after 'dummifying'
categorical variables.
}
\examples{
n <- 1000
X_df <- data.frame(income    = runif(n, 1, 1000),
                   sex       = sample(c("F", "M"), n, replace = TRUE),
                   age       = sample(100, n, replace = TRUE),
                   ethnicity = sample(LETTERS[1:5], n, replace = TRUE))

get_dummified_names(dataset = X_df)

get_dummified_names(dataset       = X_df,
                    model_formula = ~ age + ethnicity)

}
