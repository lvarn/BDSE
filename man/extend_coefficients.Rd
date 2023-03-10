% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions_synthetic_data.R
\name{extend_coefficients}
\alias{extend_coefficients}
\title{Extend a vector of coefficients}
\usage{
extend_coefficients(X, b)
}
\arguments{
\item{X}{A design matrix with column names.}

\item{b}{A named numeric vector, with names matching column names of
\code{X}.}
}
\value{
A named numeric vector of length equal to \code{ncol(X)}.
}
\description{
This function, given a design matrix \code{X} and a vector of coefficients
\code{b}, extends \code{b} with 0-valued elements to be compatible with
\code{X}, so that \code{X \%*\% b} is correctly computed.
}
\details{
This function is a helper function for \code{generate_target_dataset}. It is
used to extend the coefficients vector to match the number of columns in the
design matrix, when only some coefficients are specified (see examples).
}
\examples{
n       <- 50
dataset <- data.frame(age_grp   = factor(sample(seq(10, 65, 5), n,
                                                replace = TRUE)),
                      income    = runif(n, 1, 50000),
                      sex       = sample(c("M", "F"), n, replace = TRUE),
                      hours_grp = sample(c("[0,10)", "[10,40)", "[40,Inf)"),
                                         n, replace = TRUE))
X <- model.matrix(~ ., data = dataset)
head(X)

b <- c('(Intercept)'      = 0.5,
       age_grp15          = 2,
       age_grp35          = 1.5,
       'hours_grp[10,40)' = 3)

b_extended <- extend_coefficients(X = X, b = b)

as.matrix(b)
as.matrix(b_extended)

}
