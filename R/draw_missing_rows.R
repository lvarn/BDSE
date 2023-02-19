#' Draw from the posterior predictive distribution of the missing records
#'
#' This function draws samples from the posterior predictive distribution of the
#' covariate values of the unobserved or missing portion of the population,
#' given the data on the observed portion of the population.
#'
#' The target population dataset is \eqn{[D_{obs}, D_{mis}]}, where
#' \eqn{D_{obs}} is the observed portion (captured by at least one of the two
#' lists List 1 and List 2), and \eqn{D_{mis}} is the unobserved or missing
#' portion of the dataset (not captured by either list). The corresponding
#' covariates matrices are \eqn{X_{obs}} and \eqn{X_{mis}}, which together form
#' the covariates matrix of the target population \eqn{X = [X_{obs}, X_{mis}]}.
#'
#' To draw records (covariates vectors) for individuals not captures by the
#' lists, \code{draw_missing_data_indices} samples from the observed covariates
#' (\eqn{X_{obs}}), with probabilities adjusted using the coverage models. The
#' size of the target population, which determines the number of missing records
#' to draw, is drawn from its posterior distribution within
#' \code{draw_missing_data_indices}.
#'
#' The two design matrices \code{X_obs1} and \code{X_obs2} have the same number
#' of rows, but may have different columns (depending on the terms in their
#' corresponding logistic regression coverage models).
#'
#' @param n_observed An integer specifying the number of rows in the observed
#'   dataset.
#' @param X_obs1,X_obs2 Each a design matrix produced by applying the coverage
#'   model for the corresponding list (List 1 or List 2) to the observed dataset
#'   (\eqn{D_{obs}}).
#' @param coefficients1,coefficients2 Each a vector of individual-level effects
#'   for the logistic regression coverage model for the corresponding list (List
#'   1 or List 2).
#' @param theta A vector of probabilities (summing to 1) for the covariate
#'   distribution.
#' @param draw_n_target_posterior A function to make a draw from the posterior
#'   distribution of the target population size, given the size of the observed
#'   portion of the population.
#' @param n_target_posterior_args Arguments for the function
#'   \code{draw_n_target_posterior}.
#' @param offsets1,offsets2 Each a numeric vector of offsets in the coverage
#'   probbaility for the corresponding list; defaults to 0.
#'
#' @return A vector of indices of the sampled records (to build the unobserved
#'   or missing dataset).
#'
#' @examples
#' n_obs         <- 250  # number of observed data
#' n_clusters    <- 50   # number of clusters
#'
#' X_obs         <- cbind('(Intercept)' = rep(1, n_obs),
#'                        sexM          = rbinom(n_obs, 1, 0.45),
#'                        age           = runif(n_obs, 1, 85))
#'
#' clusters_obs  <- sample(n_clusters, n_obs, replace = TRUE)
#'
#' theta         <- 1 / (5 * n_obs) * (rep(1, n_obs) +
#'                                     rmultinom(1,
#'                                               size = 4 * n_obs,
#'                                               prob = rep(1 / n_obs, n_obs)))
#'
#' cluster_effects1 <- rnorm(n_clusters, 0, 0.5)
#' cluster_effects2 <- rnorm(n_clusters, 0, 0.5)
#'
#' offsets1 <- cluster_effects1[clusters_obs]
#' offsets2 <- cluster_effects2[clusters_obs]
#'
#' sum(theta)   # must sum to 1
#'
#' X_mis_indices <-
#'   draw_missing_data_indices(n_observed         = n_obs,
#'                             X_obs1             = X_obs,
#'                             X_obs2             = X_obs,
#'                             coefficients1      = c(0.5, -0.5, 0.01),
#'                             coefficients2      = c(0.5, 0, 0.02),
#'                             theta              = theta,
#'                             draw_n_target_posterior = draw_n_target_jeffreys,
#'                             n_target_posterior_args = NULL,
#'                             offsets1   = offsets1,
#'                             offsets2   = offsets2)
#'
#' length(X_mis_indices) + n_obs  # size of target population
#'
#' X_mis <- X_obs[X_mis_indices, ]
#' head(X_mis)
#'
#' @export
draw_missing_data_indices <- function(n_observed,
                                      X_obs1,
                                      X_obs2,
                                      coefficients1,
                                      coefficients2,
                                      theta,
                                      draw_n_target_posterior,
                                      n_target_posterior_args = NULL,
                                      offsets1                = 0,
                                      offsets2                = 0) {
  if (nrow(X_obs1) != nrow(X_obs2))
    stop("Design matrices for coverage models of Lists 1 and 2 must ",
         "have the same number of rows")

  prob_Y00_given_x <-
    calc_missing_probabilities(list(X            = X_obs1,
                                    coefficients = coefficients1,
                                    offsets      = offsets1),
                               list(X            = X_obs2,
                                    coefficients = coefficients2,
                                    offsets      = offsets2))
  
  prob_Y00_joint_x  <- prob_Y00_given_x * theta
  prob_Y00_marginal <- sum(prob_Y00_joint_x)

  ## Draw 'n_target' (size of target population);
  ## do_call() is a wrapper for do.call() that allows passing unused arguments
  n_target <- do_call(draw_n_target_posterior,
                      c(list(n_success = n_observed,
                             p_success = 1 - prob_Y00_marginal),
                        n_target_posterior_args))

  ## Row indices for 'missing' records sampled from the observed records
  sample(x       = nrow(X_obs1),
         size    = n_target - n_observed,
         replace = TRUE,
         prob    = prob_Y00_joint_x)
}


#' Draw missing counts
#' 
#' @export
draw_missing_counts <- function(X_obs1,
                                X_obs2,
                                coefficients1,
                                coefficients2,
                                theta,
                                n_observed,
                                n_target = NULL,
                                offsets1 = 0,
                                offsets2 = 0) {
  if (nrow(X_obs1) != nrow(X_obs2))
    stop("Design matrices for coverage models of Lists 1 and 2 must ",
         "have the same number of rows")

  p00_x <- calc_missing_probabilities(list(X            = X_obs1,
                                           coefficients = coefficients1,
                                           offsets      = offsets1),
                                      list(X            = X_obs2,
                                           coefficients = coefficients2,
                                           offsets      = offsets2))

  p00_x <- p00_x * theta
  p00   <- sum(p00_x)

  ## Draw n_target for x's
  if (is.null(n_target))
    n_target <- BDSE:::draw_n_target_jeffreys(n_success = n_observed,
                                              p_success = 1 - p00)

  as.vector(stats::rmultinom(n = 1, 
                             size = n_target - n_observed,
                             prob = p00_x))
}

