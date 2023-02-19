#' Draw from the posterior distribution of linear regression coefficients
#'
#' This function draws samples from the conditional posterior distribution of
#' the coefficients of a linear regression model, with known covariance: 
#' \deqn{ 
#'   Y \sim MVN(X^T \beta + \gamma, \sigma^2 I),
#' }
#' where \eqn{\gamma} is a known offset.
#' 
#' The prior distribution of the coefficients (\eqn{\beta}) is the multivariate
#' normal (MVN) distribution:
#' \deqn{
#'   MVN(\mu_0, \Sigma_0),
#' }
#' where \eqn{\mu_0} is the prior means vector and \eqn{\Sigma_0} is the prior
#' covariance matrix.
#' 
#' The posterior distribution of the coefficients is also MVN:
#' \deqn{
#'   MVN(\Sigma \gamma, \Sigma),
#' }
#' where \eqn{\Sigma = (X^T (\sigma^2 I)^{-1} X + \Sigma_0^{-1})^{-1}}, and 
#' \eqn{\gamma = X^T (\sigma^2 I)^{-1} y + \Sigma_0^{-1} \mu_0}.
#' 
#' @param X A numeric matrix of covariates (the design matrix), having dimension
#'   \eqn{m \times p}.
#' @param y A numeric vector of length \eqn{m} of observations of the response
#'   variable.
#' @param prior_mean A numeric vector of length \eqn{p} specifying the mean of
#'   the MVN prior distribution of the coefficients.
#' @param prior_cov A numeric matrix, having dimension \eqn{p \times p},
#'   specifying the covariance of the MVN prior distribution of the
#'   coefficients.
#' @param normal_sd A numeric scalar specifying the (known) standard deviation
#'   of the normal distribution.
#' @param n An integer specifying the number of draws form the posterior
#'   distribution.
#' @param offsets An optional numeric vector, having the same length as
#'   \code{y}, of known offsets (defaults to 0).
#'   
#' @return A numeric matrix, with \code{n} rows, of the sampled coefficients.
#' 
#' @examples
#' n_data       <- 1000 # number of observations
#' n_draws      <- 2    # number of posterior draws
#' X            <- model.matrix(~ .,
#'                              data.frame(x1 = runif(n_data),
#'                                         x2 = sample(letters[1:3], 
#'                                                     n_data, 
#'                                                     replace = TRUE)))
#' coefficients <- c(2, 1.3, 1.5, 0.5)
#' normal_sd    <- 2
#' offsets      <- rnorm(n_data)
#' y            <- rnorm(n_data, X %*% coefficients + offsets, normal_sd)
#' 
#' prior_mean   <- rep(0, length(coefficients))
#' prior_cov    <- 10 * diag(length(coefficients))
#' 
#' draws <- draw_coefficients_lm(X          = X,
#'                               y          = y,
#'                               prior_mean = prior_mean,
#'                               prior_cov  = prior_cov,
#'                               normal_sd  = normal_sd,
#'                               n          = n_draws,
#'                               offsets    = offsets)
#' apply(draws, 2, median)
#'
#' @export
draw_coefficients_lm <- function(X, 
                                 y,
                                 prior_mean,
                                 prior_cov,
                                 normal_sd,
                                 offsets = 0,
                                 n       = 1) {
  y <- y - offsets
  
  inv_prior_cov  <- solve(prior_cov)
  posterior_cov  <- solve((1 / normal_sd^2) * crossprod(X) + inv_prior_cov)
  
  posterior_mean <- posterior_cov %*% ((1 / normal_sd^2) * crossprod(X, y) + 
                                         inv_prior_cov %*% prior_mean)
  
  mvtnorm::rmvnorm(n, posterior_mean, posterior_cov)
}
