#' Initialise output for the Gibbs sampler
#'
#' This function initialises a list of arrays (vectors and matrices) to store
#' the results of the Gibbs samplers.
#'
#' @param n An integer specifying the number of simulations (to store parameter
#'   values for).
#' @param n_clusters An integer specifying the number of clusters.
#' @param q1,q2 Each an integer specifying the number of columns (or covariates)
#'   in the design matrix for the logistic regression coverage model of the
#'   corresponding list (List 1 or List 2).
#' @param initial_vals A named list of values to initialise the output arrays,
#'   typically a result of a call to \code{\link{draw_initial_vals}} (defaults
#'   to \code{NULL}).
#' @param n_regions An integer specifying the number of higher-level clusters
#'   (regions) in the model. (Used if higher-level cluster effects are included
#'   in the model; Defaults to \code{NULL}.)
#' @param q_clus1,q_clus2 Each an integer specifying the number of cluster-level
#'   covariates for the corresponding list (List 1 or List 2). (Used if
#'   higher-level cluster effects are included in the model; default to
#'   \code{NULL}.)
#'
#' @return A list of (NA-filled) arrays for storing parameter values.
#'
#' @examples
#' n          <- 4
#' q1 <- q2   <- 2
#' n_clusters <- 10
#' 
#' initialise_output_storage(n            = n,
#'                           q1           = q1,
#'                           q2           = q2,
#'                           initial_vals = list(coefficients1 = rep(0.5, q1),
#'                                               coefficients2 = rep(0.1, q2)))
#'
#' ## With low-level cluster effects
#' initialise_output_storage(n            = n,
#'                           q1           = q1,
#'                           q2           = q2,
#'                           n_clusters   = n_clusters,
#'                           initial_vals = 
#'                             list(coefficients1    = rep(0.5, q1),
#'                                  coefficients2    = rep(0.1, q2),
#'                                  cluster_effects1 = rep(0.25, n_clusters),
#'                                  cluster_effects2 = rep(0.01, n_clusters),
#'                                  clusters_sd1     = 1.3,
#'                                  clusters_sd2     = 2.4))
#' 
#' ## With low- and high-level cluster effects
#' q_clus1   <- 3
#' q_clus2   <- 3
#' n_regions <- 2
#'
#' initialise_output_storage(n            = n,
#'                           n_clusters   = n_clusters,
#'                           n_regions    = n_regions,
#'                           q1           = q1,
#'                           q2           = q2,
#'                           q_clus1      = q_clus1,
#'                           q_clus2      = q_clus2,
#'                           initial_vals = 
#'                             list(coefficients1      = rep(0.5, q1),
#'                                  coefficients2      = rep(0.1, q2),
#'                                  cluster_effects1   = rep(0.25, n_clusters),
#'                                  cluster_effects2   = rep(0.01, n_clusters),
#'                                  clusters_sd1       = 1.3,
#'                                  clusters_sd2       = 2.4,
#'                                  coefficients_clus1 = rep(0.1, q_clus1),
#'                                  coefficients_clus2 = rep(0.3, q_clus2),
#'                                  regions_intercept1 = 0.5,
#'                                  regions_intercept2 = 0.5,
#'                                  region_effects1    = rnorm(n_regions, 0.5, 0.8),
#'                                  region_effects2    = rnorm(n_regions, 0.5, 0.8),
#'                                  regions_sd1        = 0.8,
#'                                  regions_sd2        = 0.8))
#'  
#' @export
initialise_output_storage <- function(n,
                                      q1,
                                      q2,
                                      initial_vals = NULL,
                                      n_clusters   = NULL,
                                      n_regions    = NULL,
                                      q_clus1      = 0L,
                                      q_clus2      = 0L) {
  ## n_target and coverage params: base model
  dims <- list(n_target      = n + 1,
               coefficients1 = c(n + 1, q1),
               coefficients2 = c(n + 1, q2),
               update_coef1  = n + 1,
               update_coef2  = n + 1,
               log_like1     = n + 1,
               log_like2     = n + 1)
  
  ## Cluster params: if low-level cluster-effects included
  if (!is.null(n_clusters))
    dims <- c(dims,
              list(cluster_effects1 = c(n + 1, n_clusters),
                   cluster_effects2 = c(n + 1, n_clusters),
                   update_clus_eff1 = c(n + 1, n_clusters),
                   update_clus_eff2 = c(n + 1, n_clusters),
                   clusters_sd1     = n + 1,
                   clusters_sd2     = n + 1))
  
  ## Cluster params: if cluster-level covariates included (model L2b)
  ## Note: if ~ 1 model, mean = lm_coef
  if (q_clus1 > 0L)
    dims <- c(dims, 
              list(coefficients_clus1 = c(n + 1, q_clus1),
                   clusters_mean1     = if (q_clus1 > 1L) 
                                          c(n + 1, n_clusters) else n + 1))
  
  if (q_clus2 > 0L)
    dims <- c(dims, 
              list(coefficients_clus2 = c(n + 1, q_clus2),
                   clusters_mean2     = if (q_clus2 > 1L) 
                                          c(n + 1, n_clusters) else n + 1))
  
  ## Region params: if higher-level cluster effects included
  if (!is.null(n_regions))
    dims <- c(dims,
              list(regions_intercept1 = n + 1,
                   regions_intercept2 = n + 1,
                   region_effects1    = c(n + 1, n_regions),
                   region_effects2    = c(n + 1, n_regions),
                   regions_sd1        = n + 1,
                   regions_sd2        = n + 1))
  
  output <- lapply(dims, function(x) array(dim = x))

  for (name in names(initial_vals))
    assign_first(output[[name]], initial_vals[[name]])

  output
}


#' Draw initial parameter values for the Gibbs sampler
#'
#' This function draws random initial parameter values for the Gibbs samplers.
#' The inidividual-level coefficients and the group-level cluster effects are
#' drawn from normal distributions, while the standard deviations of the cluster
#' effects distributions are drawn from uniform distributions.
#'
#' The number of coeffcients drawn for the logistic regression coverage model
#' for each of the two lists is equal to the lengths of the mean vectors
#' \code{mean_coef1} and \code{mean_coef2}.
#'
#' If the ranges of the standard deviations for the distribution of the cluster
#' effects are not the same for both lists, then they are specified using
#' \code{sd1_min}, \code{sd1_max}, \code{sd2_min}, and \code{sd2_max}. These
#' default to \code{sd_min} and \code{sd_max}.
#'
#' The integers in the argument names refer to the indices of the two lists.
#'
#' @param n_clusters An integer representing the number of clusters in the data.
#' @param sd_min,sd_max The lower and upper bound of the uniform distribution
#'   from which the standard deviation of the distribution of the cluster
#'   effects is drawn.
#' @param sd1_min,sd1_max,sd2_min,sd2_max The lower and upper bounds of the
#'   uniform distributions from which the standard deviations of the
#'   distributions of the cluster effects are drawn, when they differ for the
#'   two lists.
#' @param mean_coef1,mean_coef2 Each a numeric vector specifying the mean of the
#'   normal distribution from which the coefficients for the corresponding list
#'   is drawn.
#' @param cov_coef1,cov_coef2 Each a matrix specifying the covariance of the
#'   normal distribution from which the coefficients for the corresponding list
#'   is drawn.
#'
#' @return A named list of initial parameter values to be passed to
#'   \code{\link{initialise_output_storage}} through
#'   \code{run_Gibbs_sampler_*} functions, with elements:
#'   \code{coefficients1},
#'   \code{coefficients2},
#'   \code{cluster_sd1},
#'   \code{cluster_sd2},
#'   \code{cluster_effects1}, and
#'   \code{cluster_effects2}.
#'
#' @examples
#' q1 <- 3
#' q2 <- 5
#'
#' ## Non-hierarchical:
#' draw_initial_vals(mean_coef1 = rep(0, q1),
#'                   mean_coef2 = rep(0, q2),
#'                   cov_coef1  = diag(1.2, nrow = q1),
#'                   cov_coef2  = diag(0.8, nrow = q2))
#' 
#' ## Hierarchical (2-level) model:
#' draw_initial_vals(mean_coef1 = rep(0, q1),
#'                   mean_coef2 = rep(0, q2),
#'                   cov_coef1  = diag(1.2, nrow = q1),
#'                   cov_coef2  = diag(0.8, nrow = q2),
#'                   n_clusters = 10,
#'                   sd_min     = 0.1,
#'                   sd_max     = 1)
#' 
#'
#' @export
## TODO: allow specifying cluster means (cluster-level LM coefficients), and
## also region effects, etc. Also, ammend for 3-level model.
draw_initial_vals <- function(mean_coef1,
                              mean_coef2,
                              cov_coef1,
                              cov_coef2,
                              n_clusters = NULL,
                              sd_min     = NULL,
                              sd_max     = NULL,
                              sd1_min    = sd_min,
                              sd1_max    = sd_max,
                              sd2_min    = sd_min,
                              sd2_max    = sd_max) {
  ## Draw initial values for coefficient vectors
  out <- list(coefficients1 = mvtnorm::rmvnorm(1, mean_coef1, cov_coef1),
              coefficients2 = mvtnorm::rmvnorm(1, mean_coef2, cov_coef2))
  
  if (!is.null(n_clusters)) {
    ## Zero-mean, constant-variance cluster effects
    clusters_sd1 <- stats::runif(1, sd1_min, sd1_max)
    clusters_sd2 <- stats::runif(1, sd2_min, sd2_max)
  
    out <- c(out,
             list(clusters_sd1     = clusters_sd1,
                  clusters_sd2     = clusters_sd2,
                  cluster_effects1 = stats::rnorm(n_clusters, sd = clusters_sd1),
                  cluster_effects2 = stats::rnorm(n_clusters, sd = clusters_sd2)))
  }
  
  out
}


#' Add storage space to an array
#' 
#' @param x An array.
#' @param n An integer.
#' 
#' @return \code{x} with \code{n} added 'rows' in the first dimension.
#' 
#' @examples
#' a <- array(1:24, dim = c(2, 3, 4))
#' BDSE:::add_storage(a, 2)
add_storage <- function(x, 
                        n) {
  abind::abind(x, array(NA, dim = c(n, dim(x)[-1])), along = 1)
}

