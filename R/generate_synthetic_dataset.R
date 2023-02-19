#' Draw records from a synthetic base dataset to create a target dataset
#'
#' This function generates a target (population) dataset by randomly sampling
#' from a given base dataset.
#' 
#' @param n_target An integer specifying the size of the target population.
#' @param base_dataset A data.frame of covariates.
#' @param weights A numeric vector of sampling weights, the same length as the
#'   number of rows in the base dataset. (Defaults to \code{NULL}.)
#' @param seed An integer specifying the seed for the random number generator.
#'   (Defaults to \code{NULL}.)
#'
#' @return A data.frame of covariates with \code{n_target} records.
#'
#' @examples
#' base_dataset <- expand.grid(sex       = c("M", "F"),
#'                             income    = seq(0, 50000, 10000),
#'                             age_grp   = as.character(1:5),
#'                             ethnicity = LETTERS[1:5])
#' 
#' target_dataset <- generate_target_dataset(n_target     = 10000,
#'                                           base_dataset = base_dataset,
#'                                           seed         = 234)
#' 
#' str(target_dataset)
#' 
#' @export
generate_target_dataset <- function(n_target,
                                    base_dataset,
                                    weights = NULL,
                                    seed    = NULL) {
  ## Check base data.frame
  check_dataframe(dataset      = base_dataset,
                  dataset_name = "base_dataset")

  if(identical(nrow(base_dataset), 0L))
    stop("'base_dataset' has 0 rows")
 
  ## Sample from base dataset to create population dataset of size 'n_target'
  if (!is.null(seed)) set.seed(seed)
 
  target_dataset <- base_dataset[sample(nrow(base_dataset), 
                                        size    = n_target,
                                        prob    = weights,
                                        replace = TRUE), ]
  
  `rownames<-`(target_dataset, seq_len(n_target))
}


#' Draw records from a synthetic base dataset to create a larger dataset
#'
#' @param n An integer specifying the size of the dataset to generate.
#' @inheritParams generate_target_dataset
#' 
#' @examples
#' base_dataset <- expand.grid(sex       = c("M", "F"),
#'                             income    = seq(0, 50000, 10000),
#'                             age_grp   = as.character(1:5),
#'                             ethnicity = LETTERS[1:5])
#' 
#' synthetic_dataset <- generate_dataset(n            = 1000,
#'                                       base_dataset = base_dataset,
#'                                       seed         = 234)
#' 
#' str(synthetic_dataset)
#' 
#' @export
generate_dataset <- function(n,
                             base_dataset,
                             weights = NULL,
                             seed    = NULL) {
  generate_target_dataset(n_target     = n,
                          base_dataset = base_dataset,
                          weights      = weights,
                          seed         = seed)
}


#' Function to apply overcoverage
#' 
#' Given a dataset, this function generates a target inclusion indicator to 
#' simulate overcoverage with respect to the target population.
#'
#' The names of the elements of the coefficients vector (\code{coefficients}) 
#' must be a subset of the column names of the design matrix created using the 
#' corresponding model formula (\code{model_formula}).
#' 
#' If a clustering variable is specified, then a \code{"cluster_target"} index
#' column is added to the dataset.
#' 
#' @inheritParams add_undercoverage
#' @param dataset A data.frame of covariates, each row representing an
#'   individual.
#' @param formula_target A formula for the fixed effects in the linear component
#'   of the logistic model for the target inclusion indicator.
#' @param coefficients_target A named numeric vector specifying the coefficients
#'   of the target inclusion logistic model.
#' @param clusters_sd_target A numeric (positive) scalar specifying the standard
#'   deviation of the cluster effects of the logistic target inclusion model. 
#'   (Defaults to \code{NULL}.)
#' @param seed An integer specifying the seed for the random number generator.
#'   (Defaults to \code{NULL}.)
#'
#' @return The input data.frame of covariates with a logical target inclusion
#'   indicator column named \code{in_target}.
#'
#' @examples
#' base_dataset   <- expand.grid(sex       = c("M", "F"),
#'                               income    = seq(0, 50000, 10000),
#'                               age_grp   = as.character(1:5),
#'                               region    = letters[1:8],
#'                               ethnicity = LETTERS[1:5])
#' 
#' union_dataset  <- generate_dataset(n            = 1000,
#'                                    base_dataset = base_dataset)
#' 
#' formula_target <- ~ sex + income + age_grp + ethnicity
#' 
#' get_dummified_names(dataset       = base_dataset, 
#'                     model_formula = formula_target)
#' 
#' coefficients_target <- c(`(Intercept)` = 0.75,
#'                          sexF          = 0.4, 
#'                          age_grp2      = 0.5,
#'                          age_grp3      = 0.3,
#'                          ethnicityC    = 0.6,
#'                          ethnicityD    = 0.5)
#' 
#' union_dataset <- add_overcoverage(dataset             = union_dataset,
#'                                   formula_target      = formula_target,
#'                                   coefficients_target = coefficients_target,
#'                                   clustering_var      = "region",
#'                                   clusters_sd_target  = 0.25)
#' head(union_dataset)
#' table(union_dataset[["in_target"]])
#' 
#' @export
add_overcoverage <- function(dataset,
                             formula_target,
                             coefficients_target,
                             clustering_var     = NULL,
                             clusters_sd_target = NULL,
                             seed               = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(dataset)
  
  ## Generate offsets for 2-level model
  if (!is.null(clustering_var)) {
    dataset[["cluster_target"]] <- match(dataset[[clustering_var]],
                                         sort(unique(dataset[[clustering_var]])))

    n_clusters      <- max(dataset[["cluster_target"]]) # integers 1:n_clusters
    cluster_effects <- stats::rnorm(n_clusters, 0, clusters_sd_target)
    offsets         <- cluster_effects[dataset[["cluster_target"]]]
  } else
      offsets <- rep(0, n)

  ## Generate target inclusion indicator
  X <- stats::model.matrix(formula_target, data = dataset)
  b <- extend_coefficients(X = X, b = coefficients_target)
  p <- invlogit(as.vector(X %*% b + offsets))
  
  dataset[["in_target"]] <- as.logical(stats::rbinom(n, size = 1, prob = p))
  
  ## Include true coverage probabilities and offsets for testing
  dataset[["offset_target"]] <- offsets
  dataset[["p_target"]]      <- p
  
  attr(dataset, "clustering_var_overcov") <- clustering_var
  
  dataset
}


#' Function to apply list undercoverage to target dataset
#'
#' This function, given coverage models for two lists, generates list inclusion
#' indicators to simulate (under)coverage in the given target (population)
#' dataset.
#'
#' The function assumes that if overcoverage is present, then a logical
#' indicator column \code{in_target} is included in \code{dataset}. This can be
#' a result of a call to \code{\link{add_overcoverage}}.
#'
#' The names of the elements of the coefficients vectors (\code{coefficients1}
#' and \code{coefficients2}) must be a subset of the column names of the design
#' matrices created using the corresponding coverage model formulas
#' (\code{formula1} and \code{formula2}).
#'
#' If a clustering variable is specified, then a \code{"cluster"} index column
#' is added to the target dataset. Currently, the same clustering variable is
#' used for both lists (when they are 2-level models).
#'
#' If overcoverage is present, it is applied to List 2, and the coverage models
#' are conditional on being included in the target (i.e. \code{in_target =
#' TRUE}).
#' 
#' If there is partial coverage (exclusions) by design (i.e. some records or
#' subpopulations are excluded from a list by design), it is applied to List 1.
#' In the presence of partial coverage, the coverage models are conditional on
#' not being excluded by design (i.e. \code{in_list = TRUE}).
#'
#' @inheritParams arguments
#' @param clustering_var A string representing the name of the clustering
#'   variable in the base dataset. (Defaults to \code{NULL}.)
#' @param coefficients1,coefficients2 Each a named numeric vector specifying the
#'   coefficients of the coverage model for the corresponding list (List 1 or
#'   List 2).
#' @param seed An integer specifying the seed for the random number generator.
#'   (Defaults to \code{NULL}.)
#' @param has_overcoverage A Boolean indicating whether the input dataset
#'   includes records of individuals not in the target dataset (defaults to
#'   \code{FALSE}).
#' @param has_partial_coverage A Boolean indicating whether the input dataset
#'   includes an indicator for partial coverage by design (defaults to
#'   \code{FALSE}).
#'
#' @return A data.frame of covariates (the input target dataset) and list
#'   inclusion indicators for the records.
#'
#' @examples
#' n_target <- 10000
#' n_union  <- n_target + 1000
#'
#' base_dataset <- expand.grid(sex       = c("M", "F"),
#'                             income    = seq(0, 50000, 10000),
#'                             age_grp   = as.character(1:5),
#'                             region    = letters[1:8],
#'                             ethnicity = LETTERS[1:5])
#'
#' ## Example 1 --------------------
#' 
#' target_dataset <- generate_target_dataset(n_target     = n_target,
#'                                           base_dataset = base_dataset,
#'                                           seed         = 200)
#'
#' union_dataset  <- rbind(target_dataset,
#'                         generate_dataset(n_union - n_target, base_dataset))
#'
#' formula_cover  <- ~ sex + income + age_grp + ethnicity
#' get_dummified_names(dataset       = base_dataset,
#'                     model_formula = formula_cover)
#'
#' coefficients   <- c(`(Intercept)` = 0.75,
#'                     sexF          = 0.5,
#'                     age_grp2      = -0.4,
#'                     age_grp3      = -0.2,
#'                     ethnicityC    = 0.5,
#'                     ethnicityD    = -0.5)
#'
#'
#' ## Case 1.1: Without overcoverage --------
#'
#' target_dataset <- add_undercoverage(dataset        = target_dataset,
#'                                     formula1       = formula_cover,
#'                                     formula2       = formula_cover,
#'                                     coefficients1  = coefficients,
#'                                     coefficients2  = coefficients,
#'                                     clustering_var = "region",
#'                                     clusters_sd1   = 0.8,
#'                                     clusters_sd2   = 0.9,
#'                                     seed           = 100)
#'
#' make_coverage_summary(target_dataset)
#'
#'
#' ## Case 1.2: With overcoverage --------
#' 
#' union_dataset[["in_target"]] <- c(rep(TRUE, n_target),
#'                                   rep(FALSE, n_union - n_target))
#'
#' union_dataset <- add_undercoverage(dataset          = union_dataset,
#'                                    formula1         = formula_cover,
#'                                    formula2         = formula_cover,
#'                                    coefficients1    = coefficients,
#'                                    coefficients2    = coefficients,
#'                                    clustering_var   = "region",
#'                                    clusters_sd1     = 0.8,
#'                                    clusters_sd2     = 0.9,
#'                                    has_overcoverage = TRUE,
#'                                    seed             = 100)
#'
#' make_coverage_summary(union_dataset)
#'
#' ## Checks
#' all(xtabs(~ L1 + L2, data = target_dataset) ==
#'     xtabs(~ L1 + L2, data = subset(union_dataset, in_target)))
#'
#' sum(is.na(union_dataset[["cluster"]])) == n_union - n_target
#'
#'
#' ## Example 2 --------------------
#' 
#' target_dataset <- generate_target_dataset(n_target     = n_target,
#'                                           base_dataset = base_dataset,
#'                                           seed         = 200)
#' union_dataset  <- rbind(target_dataset,
#'                         generate_dataset(n_union - n_target, base_dataset))
#' 
#' ## Case 2.1: With partial coverage by design and no overcoverage --------
#'
#' target_dataset[["in_list"]] <- sample(c(TRUE, FALSE), 
#'                                       n_target, 
#'                                       replace = TRUE, 
#'                                       prob    = c(0.4, 0.6))
#'                                       
#' target_dataset <- add_undercoverage(dataset        = target_dataset,
#'                                     formula1       = formula_cover,
#'                                     formula2       = formula_cover,
#'                                     coefficients1  = coefficients,
#'                                     coefficients2  = coefficients,
#'                                     clustering_var = "region",
#'                                     clusters_sd1   = 0.8,
#'                                     clusters_sd2   = 0.9,
#'                                     seed           = 100,
#'                                     has_partial_coverage = TRUE)
#'
#' make_coverage_summary(target_dataset)
#' 
#' xtabs(~ L1 + L2 + in_list, 
#'       data = data.frame(L1 = factor(target_dataset[["L1"]], 
#'                                     levels = c(TRUE, FALSE), 
#'                                     labels = c(1, 0)),
#'                         L2 = factor(target_dataset[["L2"]], 
#'                                     levels = c(TRUE, FALSE), 
#'                                     labels = c(1, 0)),
#'                         in_list = factor(target_dataset[["in_list"]], 
#'                                          levels = c(TRUE, FALSE), 
#'                                          labels = c(1, 0))))
#' 
#' ## Case 2.2: With undercoverage by design and overcoverage --------
#' 
#' union_dataset[["in_target"]] <- c(rep(TRUE, n_target),
#'                                   rep(FALSE, n_union - n_target))
#' 
#' union_dataset[["in_list"]] <- FALSE 
#' union_dataset[["in_list"]][union_dataset[["in_target"]]] <- 
#'   sample(c(TRUE, FALSE), 
#'          n_target, 
#'          replace = TRUE, 
#'          prob    = c(0.4, 0.6))
#' 
#' union_dataset <- add_undercoverage(dataset          = union_dataset,
#'                                    formula1         = formula_cover,
#'                                    formula2         = formula_cover,
#'                                    coefficients1    = coefficients,
#'                                    coefficients2    = coefficients,
#'                                    clustering_var   = "region",
#'                                    clusters_sd1     = 0.8,
#'                                    clusters_sd2     = 0.9,
#'                                    seed             = 100,
#'                                    has_overcoverage = TRUE,
#'                                    has_partial_coverage = TRUE)
#' 
#' make_coverage_summary(union_dataset)
#' 
#' xtabs(~ in_list + L1, union_dataset)
#'
#' @export
add_undercoverage <- function(dataset,
                              formula1,
                              formula2,
                              coefficients1,
                              coefficients2,
                              clustering_var   = NULL,
                              clusters_sd1     = NULL, 
                              clusters_sd2     = NULL,
                              seed             = NULL,
                              has_overcoverage = FALSE,
                              has_partial_coverage = FALSE) {
  ## Some quick checks
  check_in_target(has_overcoverage, dataset)
  check_in_list(has_partial_coverage, dataset)
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Identify List 2 records in target dataset and List 1 records to exclude
  in_target      <- if (has_overcoverage) 
                      dataset[["in_target"]] else 
                        rep(TRUE, nrow(dataset))
  
  target_dataset <- dataset[in_target, ]
  n_target       <- sum(in_target)
  
  ## Make cluster ID variable
  if (!is.null(clustering_var)) {
    target_dataset[["cluster"]] <- match(target_dataset[[clustering_var]],
                                         sort(unique(target_dataset[[clustering_var]])))

    n_clusters <- max(target_dataset[["cluster"]]) # integers 1:n_clusters
  }
  
  ## Make cluster effects and offsets for 2-level model
  cluster_effects <- lapply(c(clusters_sd1, clusters_sd2),
                            function(s) {
                              if (!is.null(s)) 
                                stats::rnorm(n_clusters, 0, s) else
                                  rep(0, n_clusters)
                            })
  
  offsets <- lapply(cluster_effects,
                    function(eff) eff[target_dataset[["cluster"]]])

  ## Generate list inclusion indicators
  P <- lapply(1:2,
              function(j) {
                X <- stats::model.matrix(get(paste0("formula", j)), 
                                         data = target_dataset)
                b <- extend_coefficients(X = X, 
                                         b = get(paste0("coefficients", j)))
                calc_probs(X, b, offsets[[j]])
              })
  
  L <- lapply(P, 
              function(p) as.logical(stats::rbinom(n_target, size = 1, prob = p)))
  
  ## Exclude target records not in List 1 by design 
  if (has_partial_coverage) {
    P[[1]][!target_dataset[["in_list"]]] <- NA
    L[[1]][!target_dataset[["in_list"]]] <- FALSE
  }
  
  ## Extend target dataset to include list indicators
  target_dataset[["L1"]]  <- L[[1]]
  target_dataset[["L2"]]  <- L[[2]]
  target_dataset[["Y11"]] <-  target_dataset[["L1"]] &  target_dataset[["L2"]]
  target_dataset[["Y10"]] <-  target_dataset[["L1"]] & !target_dataset[["L2"]]
  target_dataset[["Y01"]] <- !target_dataset[["L1"]] &  target_dataset[["L2"]]
  target_dataset[["Y00"]] <- !target_dataset[["L1"]] & !target_dataset[["L2"]]
  
  ## Include true coverage probabilities and offsets for testing
  target_dataset[["offset1"]]  <- offsets[[1]]
  target_dataset[["offset2"]]  <- offsets[[2]]
  target_dataset[["p1"]]       <- P[[1]]
  target_dataset[["p2"]]       <- P[[2]]
  
  ## Assign NA to undercoverage fields for List 2 records not in target
  if (has_overcoverage) {
    dataset[setdiff(names(target_dataset), names(dataset))] <- NA
    dataset[in_target, ] <- target_dataset

    ## Update cell/list indicators for records outside target
    indicators <- grep("^L\\d{1}$|^Y\\d{2}$", names(dataset), value = TRUE)
    dataset[!in_target, indicators] <- FALSE
    dataset[!in_target, "L2"]       <- TRUE  # outside target (overcoverage)
  } else
      dataset <- target_dataset
  
  attr(dataset, "clustering_var") <- clustering_var

  dataset
}


## TODO: incorporate partial coverage in LE functions.

#' Function to apply linkage error to observed dataset
#' 
#' This function generates missed link indicators to simulate linkage error in
#' the given observed dataset: the observed dataset is extended with the matches
#' for the missed links.
#' 
#' Specificity is assumed to be 1 (i.e. no false positive links).
#' 
#' If \code{n_validated} is specified, then linkage error correction is made to
#' a random sample of size \code{n_validated}, to simulate the linkage
#' validation study.
#'
#' The names of the elements of the coefficients vector (\code{coefficients})
#' must be a subset of the column names of the design matrix created using the
#' sensitivity model formula (\code{formula}).
#' 
#' @param dataset A data.frame of covariates and cell indicators. (This is the
#'   observed data, i.e. it does not include records from cell (0,0).)
#' @param formula_sens A formula, the sensitivity model (used for generating the
#'   observed links binary indicator).
#' @param coefficients_sens A named numeric vector specifying the coefficients
#'   of the sensitivity model.
#' @param seed An integer specifying the seed for the random number generator.
#'   (Defaults to \code{NULL}.)
#' @param n_validated An integer specifying the size of the validation sample.
#'   (Defaults to \code{NULL}.)
#' 
#' @return A data.frame of covariates and cell indicators with linkage error.
#'   (True cell and list indicators are also returned only for the purpose of
#'   comparison/testing.)
#'
#' @examples
#' n_target <- 10000
#' n_union  <- n_target + 1000
#' 
#' base_dataset <- expand.grid(sex       = c("M", "F"),
#'                             income    = seq(0, 50000, 10000),
#'                             age_grp   = as.character(1:5),
#'                             region    = letters[1:8],
#'                             ethnicity = LETTERS[1:5])
#' 
#' union_dataset <- generate_dataset(n            = n_union,
#'                                   base_dataset = base_dataset)
#' 
#' union_dataset[["in_target"]] <- c(rep(TRUE, n_target), 
#'                                   rep(FALSE, n_union - n_target))
#'
#' formula_cover  <- ~ sex + income + age_grp + ethnicity
#' get_dummified_names(dataset       = base_dataset, 
#'                     model_formula = formula_cover)
#' 
#' coefficients   <- c(`(Intercept)` = 0.1,
#'                     sexF          = 0.75, 
#'                     age_grp2      = -0.4,
#'                     age_grp3      = -0.2,
#'                     ethnicityC    = 1.2,
#'                     ethnicityD    = -0.5)
#' 
#' union_dataset <- add_undercoverage(dataset           = union_dataset,
#'                                    formula1          = formula_cover,
#'                                    formula2          = formula_cover, 
#'                                    coefficients1     = coefficients,
#'                                    coefficients2     = coefficients,
#'                                    clustering_var    = "region",
#'                                    clusters_sd1      = 0.8,
#'                                    clusters_sd2      = 0.9,
#'                                    has_overcoverage  = TRUE)
#' 
#' make_coverage_summary(union_dataset)
#' 
#' formula_sens      <- ~ sex + ethnicity
#' coefficients_sens <- c(`(Intercept)` = 0.8,
#'                        sexF          = 0.1)
#' 
#' obs_dataset    <- union_dataset[!union_dataset[["Y00"]], ]
#' obs_dataset_LE <- add_linkage_error(dataset           = obs_dataset,
#'                                     formula_sens      = formula_sens,
#'                                     coefficients_sens = coefficients_sens,
#'                                     n_validated       = 1000)
#'                                     
#' ## Checks
#' xtabs(~ Y11, data = subset(obs_dataset, Y11 | Y10))
#' xtabs(~ YLE11 + Y11, data = subset(obs_dataset_LE, YLE11 | YLE10))
#' 
#' all(xtabs(~ Y11, data = subset(obs_dataset, Y11 | Y10)) ==
#'     colSums(xtabs(~ YLE11 + Y11, data = subset(obs_dataset_LE, YLE11 | YLE10))))
#' 
#' @export
add_linkage_error <- function(dataset,
                              formula_sens,
                              coefficients_sens,
                              seed              = NULL,
                              n_validated       = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n_obs <- nrow(dataset)  # before adding LE

  ## Sensitivity: p = P{YLE11 = 1 | Y11 = 1, X}
  true_links <- which(dataset[["Y11"]])
 
  X <- stats::model.matrix(formula_sens, data = dataset[true_links, ])
  b <- extend_coefficients(X = X, b = coefficients_sens)
  p <- invlogit(as.vector(X %*% b))

  ## Initialize YLEs with true cell indicator values
  dataset[["YLE11"]] <- dataset[["Y11"]]
  dataset[["YLE10"]] <- dataset[["Y10"]]
  dataset[["YLE01"]] <- dataset[["Y01"]]
  dataset[["YLE00"]] <- dataset[["Y00"]]
  
  ## Generate missed links
  dataset[["YLE11"]][true_links] <- 
    as.logical(stats::rbinom(length(true_links), size = 1, prob = p))

  missed_links <- which(dataset[["Y11"]] & !dataset[["YLE11"]])
  dataset[["YLE10"]][missed_links] <- TRUE
  
  ## Simulate validation study
  dataset[["validated"]] <- FALSE

  if (!is.null(n_validated)) {
    if (n_validated > n_obs)
      stop("n_validated must be smaller than nrow(dataset)")
    
    list1 <- which(dataset[["Y11"]] | dataset[["Y10"]])
    dataset[["validated"]][sample(list1, n_validated)] <- TRUE
  }

  ## Duplicate, in the (0,1) cell, all (1,0) records marked as missed links,
  ## excluding validated records
  missed_matches <- intersect(missed_links, which(!dataset[["validated"]]))
  dataset        <- rbind(dataset, dataset[missed_matches, ])
  
  dataset[["YLE10"]][(n_obs + 1):nrow(dataset)] <- FALSE
  dataset[["YLE01"]][(n_obs + 1):nrow(dataset)] <- TRUE

  dataset
}


#' Make observed dataset
#'
#' Makes the observed dataset for the Gibbs sampler by excluding records not on
#' either of the two lists and creating a sequential clusters column if some
#' clusters in the target dataset are not observed.
#'
#' @param target_dataset A data.frame, the target dataset with cell indicator
#'   and cluster columns.
#' @param clustering_var A string specifying the clustering variable (the name
#'   of a column in target dataset).
#' @param clusters_col A string specifying the sequential cluster IDs column in
#'   the target dataset. (Defaults to "cluster".)
#'
#' @return A data.frame, the observed dataset. If some clusters are not
#'   observed, then an additional sequential cluster ID column called
#'   \code{paste0(clusters_col, "_obs")} will be added.
#' 
#' @export
make_observed_dataset <- function(target_dataset,
                                  clustering_var,
                                  clusters_col = "cluster") {
  obs_dataset <- target_dataset[!target_dataset[["Y00"]], ]
  
  ## Check observed clusters and make a sequential variant
  clusters_target <- unique(target_dataset[[clustering_var]])
  clusters_obs    <- unique(obs_dataset[[clustering_var]])
  
  if ( !all(clusters_target %in% clusters_obs) ) {
    obs_dataset[[paste0(clusters_col, "_obs")]] <- 
      match(obs_dataset[[clusters_col]], 
            sort(unique(obs_dataset[[clusters_col]])))
  }
  
  attr(obs_dataset, "clusters_not_obs") <- 
    setdiff(clusters_target, clusters_obs)
  
  obs_dataset
}


#' Make cluster-level dataset
#'
#' Makes the cluster-level (area-level) dataset indicating observed and
#' unobserved clusters.
#'
#' @param target_dataset,obs_dataset Each s data.frame, the target or observed
#'   dataset including the clustering variables.
#' @param joining_vars A character vector specifying the clustering variables on
#'   which to join the target and observed datasets.
#' @param target_vars,obs_vars Each a character vector specifying any additional
#'   target or observed dataset columns to be included in the output. (Both
#'   default to \code{NULL}.)
#' @param target_clusters_col,obs_clusters_col Each a string specifying the
#'   sequential cluster IDs column in the target or observed dataset. (Defaults
#'   to "cluster" for the target dataset and "cluster_obs" for the observed
#'   dataset.)
#'
#' @return A data.frame, the merged target and observed cluster-level datasets.
#'   The observed variable columns will be \code{NA} for all clusters in the
#'   target dataset that were not observed.
#'
#' @export
## TODO: test this and add bit for splines!
make_clusters_dataset <- function(target_dataset,
                                  obs_dataset,
                                  joining_vars,
                                  target_vars = NULL,
                                  obs_vars    = NULL,
                                  target_clusters_col = "cluster",
                                  obs_clusters_col    = "cluster_obs") {
  ## Target cluster-level dataset with cluster frequencies
  clus_cols <- c(joining_vars, target_clusters_col, target_vars)
 
  target_clus_dataset <- distribute(cbind(target_dataset[clus_cols], count = 1),
                                    contract_cols = "count")

  ## Observed cluster-level dataset
  clus_rows_obs <- match(seq_along(unique(obs_dataset[[obs_clusters_col]])), 
                         obs_dataset[[obs_clusters_col]])
  clus_cols_obs <- c(joining_vars, obs_clusters_col, obs_vars)
  
  obs_clus_dataset <- obs_dataset[clus_rows_obs, clus_cols_obs]
  
  ## Merge with observed cluster-level dataset
  clus_dataset <- merge(x = target_clus_dataset, 
                        y = obs_clus_dataset,
                        by = joining_vars,
                        all.x = TRUE,
                        suffixes = c("", "_obs"))
  
  clus_dataset
}



