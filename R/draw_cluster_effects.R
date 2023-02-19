#' Draw from the posterior distribution of group-level effects
#'
#' This function uses the Metropolis-Hastings algorithm with the symmetric
#' proposal distribution
#' \preformatted{
#'  Normal(current_cluster_effects, proposal_sd),
#' }
#' to draw from the posterior distribution of the group-level or cluster-level
#' effects. These parameters form the offsets (or random intercept) of a
#' logistic regression model. The prior distribution of these cluster effects,
#' with known mean and standard deviation, is given by
#' \preformatted{
#'  Normal(prior_mean, prior_sd).
#' }
#'
#' This function is a single iteration of the Metropolis-Hastings algorithm,
#' and is intended for use within a Gibbs sampler.
#'
#' @inheritParams arguments
#' @param current_cluster_effects A numeric vector giving current cluster effect
#'   values.
#' @param prior_mean A numeric vector giving the means for the Normal prior
#'   distribution on the cluster effects.
#' @param prior_sd The standard deviation for the Normal prior distribution on
#'   the cluster effects.
#' @param proposal_sd A numeric value giving the standard deviation for the
#'   (Normal) proposal distribution, with mean vector
#'   \code{current_cluster_effects}.
#' @param offsets A numeric vector giving the known offsets at the unit record
#'   level.
#'
#' @return A named list, with the following elements:
#' \itemize{
#'     \item \code{cluster_effects}: the vector of sampled cluster_effects;
#'     \item \code{was_updated}: a logical vector, the same length as
#'       \code{cluster_effects}, indicating whether the cluster effects
#'       were updated in the Metropolis-Hastings run.
#' }
#'
#' @examples
#' n               <- 1000    # number of observations
#' n_clusters      <- 50      # number of clusters
#' clusters        <- sample(n_clusters, n, replace = TRUE)
#' X               <- cbind('(Intercept)' = rep(1, n),
#'                          age           = runif(n, 1, 85),
#'                          sexM          = rbinom(n, 1, 0.45))
#' y               <- rbinom(n, 1, 0.7)
#' coefficients    <- c(0.5, 0.01, -1)         # individual-level effects
#' cluster_effects <- runif(n_clusters, -1, 1) # group-level effects
#' prior_mean      <- rnorm(n_clusters, 0, 1)
#' prior_sd        <- 2
#' proposal_sd     <- 1.25
#' offsets         <- rnorm(n, 0, 1)  # unit-record-level offsets
#'
#' res <- draw_cluster_effects(current_cluster_effects = cluster_effects,
#'                             prior_mean              = prior_mean,
#'                             prior_sd                = prior_sd,
#'                             proposal_sd             = proposal_sd,
#'                             X                       = X,
#'                             y                       = y,
#'                             coefficients            = coefficients,
#'                             clusters                = clusters,
#'                             n_clusters              = n_clusters,
#'                             offsets                 = offsets)
#'
#' head(res[["cluster_effects"]])
#' table(res[["was_updated"]])
#'
#' @export
draw_cluster_effects <- function(current_cluster_effects,
                                 prior_mean,
                                 prior_sd,
                                 proposal_sd,
                                 X,
                                 y,
                                 coefficients,
                                 clusters,
                                 n_clusters,
                                 offsets = 0,
                                 weights = 1,
                                 counts  = 1,
                                 min = -Inf,
                                 max = Inf) {
  if(1L < length(offsets) & length(offsets) < nrow(X))
    stop("'offsets' must be a scalar or a vector of length nrow(X)")
 
  cluster_effects_alternatives <- 
    list(current  = current_cluster_effects,
         proposed = rtruncated_norm(n_clusters, 
                                    current_cluster_effects,
                                    proposal_sd,
                                    min = min, 
                                    max = max))

  log_prior <- lapply(cluster_effects_alternatives,
                      log_dnorm,
                      mean = prior_mean,
                      sd   = prior_sd)

  offsets_alternatives <- lapply(cluster_effects_alternatives,
                                 function(cluster_effects)
                                   cluster_effects[clusters] + offsets)

  log_likelihood <- lapply(offsets_alternatives,
                           function(offsets)
                             calc_loglike_logistic(X            = X,
                                                   y            = y,
                                                   coefficients = coefficients,
                                                   clusters     = clusters,
                                                   n_clusters   = n_clusters,
                                                   offsets      = offsets,
                                                   weights      = weights,
                                                   counts       = counts))

  ## Since the proposal distribution is symmetric, the log MH-ratio is
  log_ratio <- log_likelihood[["proposed"]] + log_prior[["proposed"]] -
    (log_likelihood[["current"]] + log_prior[["current"]])

  accept <- is.nan(log_ratio) | log(stats::runif(n_clusters)) < log_ratio

  list(cluster_effects = ifelse(accept,
                                cluster_effects_alternatives[["proposed"]],
                                cluster_effects_alternatives[["current"]]),
       was_updated     = accept)
}

