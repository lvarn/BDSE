% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_synthetic_dataset.R
\name{add_undercoverage}
\alias{add_undercoverage}
\title{Function to apply list undercoverage to target dataset}
\usage{
add_undercoverage(
  dataset,
  formula1,
  formula2,
  coefficients1,
  coefficients2,
  clustering_var = NULL,
  clusters_sd1 = NULL,
  clusters_sd2 = NULL,
  seed = NULL,
  has_overcoverage = FALSE,
  has_partial_coverage = FALSE
)
}
\arguments{
\item{formula1, formula2}{Each a formula for creating the design matrix for
the coverage model of the corresponding list (list 1 or 2).}

\item{coefficients1, coefficients2}{Each a named numeric vector specifying the
coefficients of the coverage model for the corresponding list (List 1 or
List 2).}

\item{clustering_var}{A string representing the name of the clustering
variable in the base dataset. (Defaults to \code{NULL}.)}

\item{seed}{An integer specifying the seed for the random number generator.
(Defaults to \code{NULL}.)}

\item{has_overcoverage}{A Boolean indicating whether the input dataset
includes records of individuals not in the target dataset (defaults to
\code{FALSE}).}

\item{has_partial_coverage}{A Boolean indicating whether the input dataset
includes an indicator for partial coverage by design (defaults to
\code{FALSE}).}
}
\value{
A data.frame of covariates (the input target dataset) and list
  inclusion indicators for the records.
}
\description{
This function, given coverage models for two lists, generates list inclusion
indicators to simulate (under)coverage in the given target (population)
dataset.
}
\details{
The function assumes that if overcoverage is present, then a logical
indicator column \code{in_target} is included in \code{dataset}. This can be
a result of a call to \code{\link{add_overcoverage}}.

The names of the elements of the coefficients vectors (\code{coefficients1}
and \code{coefficients2}) must be a subset of the column names of the design
matrices created using the corresponding coverage model formulas
(\code{formula1} and \code{formula2}).

If a clustering variable is specified, then a \code{"cluster"} index column
is added to the target dataset. Currently, the same clustering variable is
used for both lists (when they are 2-level models).

If overcoverage is present, it is applied to List 2, and the coverage models
are conditional on being included in the target (i.e. \code{in_target =
TRUE}).

If there is partial coverage (exclusions) by design (i.e. some records or
subpopulations are excluded from a list by design), it is applied to List 1.
In the presence of partial coverage, the coverage models are conditional on
not being excluded by design (i.e. \code{in_list = TRUE}).
}
\examples{
n_target <- 10000
n_union  <- n_target + 1000

base_dataset <- expand.grid(sex       = c("M", "F"),
                            income    = seq(0, 50000, 10000),
                            age_grp   = as.character(1:5),
                            region    = letters[1:8],
                            ethnicity = LETTERS[1:5])

## Example 1 --------------------

target_dataset <- generate_target_dataset(n_target     = n_target,
                                          base_dataset = base_dataset,
                                          seed         = 200)

union_dataset  <- rbind(target_dataset,
                        generate_dataset(n_union - n_target, base_dataset))

formula_cover  <- ~ sex + income + age_grp + ethnicity
get_dummified_names(dataset       = base_dataset,
                    model_formula = formula_cover)

coefficients   <- c(`(Intercept)` = 0.75,
                    sexF          = 0.5,
                    age_grp2      = -0.4,
                    age_grp3      = -0.2,
                    ethnicityC    = 0.5,
                    ethnicityD    = -0.5)


## Case 1.1: Without overcoverage --------

target_dataset <- add_undercoverage(dataset        = target_dataset,
                                    formula1       = formula_cover,
                                    formula2       = formula_cover,
                                    coefficients1  = coefficients,
                                    coefficients2  = coefficients,
                                    clustering_var = "region",
                                    clusters_sd1   = 0.8,
                                    clusters_sd2   = 0.9,
                                    seed           = 100)

make_coverage_summary(target_dataset)


## Case 1.2: With overcoverage --------

union_dataset[["in_target"]] <- c(rep(TRUE, n_target),
                                  rep(FALSE, n_union - n_target))

union_dataset <- add_undercoverage(dataset          = union_dataset,
                                   formula1         = formula_cover,
                                   formula2         = formula_cover,
                                   coefficients1    = coefficients,
                                   coefficients2    = coefficients,
                                   clustering_var   = "region",
                                   clusters_sd1     = 0.8,
                                   clusters_sd2     = 0.9,
                                   has_overcoverage = TRUE,
                                   seed             = 100)

make_coverage_summary(union_dataset)

## Checks
all(xtabs(~ L1 + L2, data = target_dataset) ==
    xtabs(~ L1 + L2, data = subset(union_dataset, in_target)))

sum(is.na(union_dataset[["cluster"]])) == n_union - n_target


## Example 2 --------------------

target_dataset <- generate_target_dataset(n_target     = n_target,
                                          base_dataset = base_dataset,
                                          seed         = 200)
union_dataset  <- rbind(target_dataset,
                        generate_dataset(n_union - n_target, base_dataset))

## Case 2.1: With partial coverage by design and no overcoverage --------

target_dataset[["in_list"]] <- sample(c(TRUE, FALSE), 
                                      n_target, 
                                      replace = TRUE, 
                                      prob    = c(0.4, 0.6))
                                      
target_dataset <- add_undercoverage(dataset        = target_dataset,
                                    formula1       = formula_cover,
                                    formula2       = formula_cover,
                                    coefficients1  = coefficients,
                                    coefficients2  = coefficients,
                                    clustering_var = "region",
                                    clusters_sd1   = 0.8,
                                    clusters_sd2   = 0.9,
                                    seed           = 100,
                                    has_partial_coverage = TRUE)

make_coverage_summary(target_dataset)

xtabs(~ L1 + L2 + in_list, 
      data = data.frame(L1 = factor(target_dataset[["L1"]], 
                                    levels = c(TRUE, FALSE), 
                                    labels = c(1, 0)),
                        L2 = factor(target_dataset[["L2"]], 
                                    levels = c(TRUE, FALSE), 
                                    labels = c(1, 0)),
                        in_list = factor(target_dataset[["in_list"]], 
                                         levels = c(TRUE, FALSE), 
                                         labels = c(1, 0))))

## Case 2.2: With undercoverage by design and overcoverage --------

union_dataset[["in_target"]] <- c(rep(TRUE, n_target),
                                  rep(FALSE, n_union - n_target))

union_dataset[["in_list"]] <- FALSE 
union_dataset[["in_list"]][union_dataset[["in_target"]]] <- 
  sample(c(TRUE, FALSE), 
         n_target, 
         replace = TRUE, 
         prob    = c(0.4, 0.6))

union_dataset <- add_undercoverage(dataset          = union_dataset,
                                   formula1         = formula_cover,
                                   formula2         = formula_cover,
                                   coefficients1    = coefficients,
                                   coefficients2    = coefficients,
                                   clustering_var   = "region",
                                   clusters_sd1     = 0.8,
                                   clusters_sd2     = 0.9,
                                   seed             = 100,
                                   has_overcoverage = TRUE,
                                   has_partial_coverage = TRUE)

make_coverage_summary(union_dataset)

xtabs(~ in_list + L1, union_dataset)

}
