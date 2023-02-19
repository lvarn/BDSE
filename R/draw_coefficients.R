#' Draw from the posterior distribution of individual-level effects
#'
#' This function uses the Metropolis-Hastings algorithm with the symmetric MVN
#' distribution
#' \preformatted{
#'   MVN(current_coefficients, proposal_cov),
#' }
#' to draw from the posterior distribution of the individual-level effects of
#' a logistic regression model. The prior distribution, with known mean and
#' covariance, is given by
#' \preformatted{
#'   MVN(prior_mean, prior_cov).
#' }
#'
#' This function is a single iteration of the Metropolis-Hastings algorithm,
#' and is intended for use within a Gibbs sampler.
#'
#'
#' @inheritParams arguments
#' @param current_coefficients A numeric vector giving current coefficient
#' values (individual-level effects).
#' @param prior_mean A vector giving the mean for the MVN prior distribution
#'   on the coefficients.
#' @param prior_cov A matrix giving the covariance for the MVN prior
#'   distribution on the coefficients.
#' @param proposal_cov A matrix giving the covariance for the MVN proposal
#'   distribution, with mean vector \code{current_coefficients}.
#' @param proposal_scale A numeric value to scale the covariance matrix
#'   \code{proposal_cov} of the proporsal distribution.
#'
#' @return A named list, with the following elements:
#' \itemize{
#'     \item \code{coefficients}: the vector of sampled coefficients;
#'     \item \code{log_prior}: the prior log-likelihood of the sampled
#'             coefficient values;
#'     \item \code{log_likelihood}: the log-likelihood of the data model for
#'             the sampled coefficient values; and
#'     \item \code{was_updated}: a logical indicator of whether the coefficients
#'             vector was updated in the Metropolis-Hastings run.
#' }
#'
#' @examples
#' n               <- 1000  # number of observations
#' n_clusters      <- 50    # number of clusters
#' clusters        <- sample(n_clusters, n, replace = TRUE)
#' X               <- cbind('(Intercept)' = rep(1, n),
#'                          age           = runif(n, 1, 85))
#' coeffs          <- c('(Intercept)' = 0.5, 
#'                      age           = 0.05)
#' y               <- rbinom(n, 1, 0.7)
#' cluster_effects <- rnorm(n_clusters, 0, 2)
#' offsets         <- cluster_effects[clusters]
#' 
#' draw_coefficients(current_coefficients = coeffs,
#'                   prior_mean           = rep(0, times = 2),
#'                   prior_cov            = diag(5^2, nrow = 2),
#'                   proposal_cov         = matrix(c(1.5^2, 1, 1, 1.5^2), nrow = 2),
#'                   X                    = X,
#'                   y                    = y,
#'                   offsets              = offsets)
#' 
#' @export
draw_coefficients <- function(current_coefficients,
                              prior_mean,
                              prior_cov,
                              proposal_cov,
                              X,
                              y,
                              offsets        = 0,
                              weights        = 1,
                              counts         = 1,
                              proposal_scale = 1) {
  coefficient_alternatives <-
    list(current  = current_coefficients,
         proposed = as.vector(
                      mvtnorm::rmvnorm(n      = 1,
                                       mean   = current_coefficients,
                                       sigma  = proposal_cov * proposal_scale,
                                       method = "chol")))

  log_prior <- lapply(coefficient_alternatives,
                      mvtnorm::dmvnorm,
                      mean  = prior_mean,
                      sigma = as.matrix(prior_cov),
                      log   = TRUE)

  log_likelihood <- lapply(coefficient_alternatives,
                           function(b)
                             calc_loglike_logistic(y            = y,
                                                   X            = X,
                                                   coefficients = b,
                                                   offsets      = offsets,
                                                   weights      = weights,
                                                   counts       = counts))

  log_ratio <- log_likelihood[["proposed"]] + log_prior[["proposed"]] -
    (log_likelihood[["current"]] + log_prior[["current"]])

  accept <- is.nan(log_ratio) | log(stats::runif(1)) < log_ratio
  out    <- if (accept) "proposed" else "current"

  list(coefficients      = coefficient_alternatives[[out]],
       log_prior         = log_prior[[out]],
       log_likelihood    = log_likelihood[[out]],
       was_updated       = accept)
}

