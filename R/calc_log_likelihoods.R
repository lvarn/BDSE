#' Calculate log-likelihood of binomial number of trials
#'
#' This helper function is used in \code{\link[BDSE]{draw_n_target_uniform}}, to
#' calculate the log-likelihood of the unknown number of trials \code{N} of a
#' binomial model, for known number of successes and probability of success.
#'
#' @param N An integer specifying the number of trials.
#' @param k An integer specifying the number of successes.
#' @param p A numeric value in \eqn{(0, 1)} specifying the probability of
#'   success in each trial.
#'
#' @return A numeric value (the logarithm of the likelihood).
#'
#' @examples
#' calc_loglike_binom_n_trial(N = 250,
#'                            k = 200,
#'                            p = 0.6)
#'
#' @export
calc_loglike_binom_n_trial <- function(N, k, p) {
  lfactorial(N) - (lfactorial(N - k) + lfactorial(k)) + (N - k) * log(1 - p)
}


#' Calculate log-likelihood of a logistic regression model
#'
#' This function calculates the log-likelihood for a vector of observed
#' Bernoulli indicators, for a logistic regression model:
#' \deqn{
#'    [Y | \beta_0, \beta_1, X = x] ~ Bernoulli(p(x; \beta_0, \beta_1)) ,
#' }
#' where \eqn{\beta_0} is a group-level effect and \eqn{\beta_1} is a
#' vector of individual-level effects. The logistic model is
#' \deqn{
#'   logit Pr(Y = 1 | X = x, \beta_0, \beta_1) = \beta_0 + x^T \beta_1 .
#' }
#'
#' In \code{calc_loglike_logistic}, \code{coefficients} is the vector of
#' individual-level effects \eqn{\beta_1}, and
#' \code{offsets} is the vector of group-level (or cluster-level) effects
#' \eqn{\beta_0}.
#'
#' If \code{clusters} is given, then the log-likelihood values are calculated
#' for each cluster.
#'
#' If the data does not include observations from all clusters, i.e. the
#' number of unique values in the vector of observed clusters \code{clusters}
#' is less than \code{n_clusters}, then the log-likelihood for the unobserved
#' clusters is set to zero.
#'
#' @inheritParams arguments
#'
#' @return If \code{is.null(clusters)}, a scalar log-likelihood. If clusters
#' are provided, a named vector of length \code{n_clusters} of log-likelihoods
#' by cluster.
#'
#' @examples
#' n            <- 100   # number of observations
#' n_clusters   <- 5    # number of clusters
#' clusters     <- sample(n_clusters - 1, n, replace = TRUE)
#' 
#' X            <- cbind('(Intercept)' = rep(1, n),
#'                       age  = runif(n, 1, 85),
#'                       sexM = rbinom(n, 1, 0.8))
#' coefficients <- c('(Intercept)' = 0.2,
#'                   age           = 0.01,
#'                   sexM          = 0.75)
#' y            <- rbinom(n, 1, 0.03)
#' offsets      <- rnorm(n_clusters, 0, 1)[clusters]
#' weights      <- sample(3, n, replace = TRUE)
#'
#' ll1 <- BDSE:::calc_loglike_logistic(X            = X,
#'                                     y            = y,
#'                                     coefficients = coefficients, 
#'                                     clusters     = clusters,
#'                                     n_clusters   = n_clusters,
#'                                     offsets      = offsets,
#'                                     weights      = weights)
#' ll1
#'
#' ll2 <- BDSE:::calc_loglike_logistic(X            = X,
#'                                     y            = y,
#'                                     coefficients = coefficients, 
#'                                     offsets      = offsets,
#'                                     weights      = weights)
#' 
#' ll2
#' 
#' sum(ll1) == ll2
#' 
#' @export
calc_loglike_logistic <- function(X,
                                  y,
                                  coefficients,
                                  clusters   = NULL,
                                  n_clusters = NULL,
                                  offsets    = 0,
                                  weights    = 1,
                                  counts     = 1) {
  ## Log-likelihoods
  Xb   <- X %*% coefficients + offsets
  prob <- invlogit(as.vector(Xb))
  log_likelihoods <- weights * (y * Xb + counts * log(1 - prob))

  ## If clusters is provided, we return log-likelihood values by cluster
  if (!is.null(clusters)) {
    cluster_log_likelihoods <- t_apply(log_likelihoods, clusters, sum)

    ## Set 0 log-likelihood values for unrepresented clusters
    n_clusters_obs <- length(unique(clusters))
    if (!is.null(n_clusters) & n_clusters_obs < n_clusters) {
      cluster_log_likelihoods_mis        <- rep(0, n_clusters - n_clusters_obs)
      names(cluster_log_likelihoods_mis) <- setdiff(1:n_clusters,
                                                    unique(clusters))
      cluster_log_likelihoods <- order_by_name(c(cluster_log_likelihoods,
                                                 cluster_log_likelihoods_mis))
    }
    cluster_log_likelihoods
  } else {
    sum(log_likelihoods)
  }
}

