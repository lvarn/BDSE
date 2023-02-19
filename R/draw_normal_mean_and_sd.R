#' Draw from conditional posterior of the mean of a normal data model
#' 
#' This function draws samples from the posterior distribution of the mean
#' parameter of a normal model, given the variance parameter:
#' \deqn{
#'   Y \sim Normal(\mu, \sigma^2),
#' } 
#' for the conjugate normal prior on the mean parameter \eqn{\mu}:
#' \deqn{
#'   Normal(\mu_0, \sigma_0^2).
#' }
#' 
#' @param y A numeric vector of observations (data).
#' @param normal_sd A numeric scalar specifying the standard deviation of the
#'   normal model.
#' @param prior_mean A numeric scalar specifying the prior mean of the the mean
#'   parameter of the normal distribution.
#' @param prior_sd A numeric scalar specifying the prior standard deviation of
#'   the mean parameter of the normal distribution.
#' @param n An integer specifying the number of draws from the distribution.
#'   (Defaults to 1.)
#'   
#' @return A numeric vector of length \code{n} of sampled means.
#' 
#' @examples
#' draw_mean_normal(y          = rnorm(100, 3, 0.5),
#'                  normal_sd  = 2,
#'                  prior_mean = 0,
#'                  prior_sd   = 4,
#'                  n          = 5)
#' @name draw_mean
NULL


#' @rdname draw_mean
#'
#' @export
draw_mean_normal <- function(y,
                             normal_sd, 
                             prior_mean, 
                             prior_sd, 
                             n = 1) {
  posterior_var  <- prior_sd^2 * normal_sd^2 / 
                        (length(y) * prior_sd^2 + normal_sd^2)
  posterior_mean <- posterior_var * (sum(y) / normal_sd^2 + 
                                       prior_mean / prior_sd^2)
  
  stats::rnorm(n, posterior_mean, sqrt(posterior_var))
}


#' Sample from posterior distribution of the standard deviation parameter of a 
#' Normal model
#'
#' Generate draws from the posterior for the standard deviation of the 
#' normal model with known mean(s):
#' \deqn{
#'    Y_i ~ Normal(\mu_i, \sigma^2) .
#'  }
#' This function is used to update the standard deviation parameter of the 
#' normal prior distribution of the cluster effects.
#' 
#' The draws are made from the posterior distribution of either the precision 
#' \eqn{\tau^2 = 1 / \sigma^2}, or the variance \eqn{\sigma^2}, and then 
#' transformed to get the standard deviation.
#' 
#' \code{draw_sd_uniform} returns posterior draws of \eqn{\sigma} when 
#' the prior is \eqn{\sigma ~ U(prior_lower, prior_upper)}. The draws are made 
#' from the posterior distribution of the precision, which is the Gamma 
#' distribution
#' \deqn{
#'    Gamma ( (m - 1) / 2, \sum(Y_i - \mu_i)^2 / 2 ),
#' }
#' truncated to include only \eqn{\sigma \in (prior_lower, prior_upper)}, 
#' where \eqn{m} is the number of residuals. The drawn precisions are then 
#' transformed to get the standard deviations: \eqn{\sigma = 1 / \sqrt \tau^2}.
#' 
#' \code{draw_sd_invchisq} returns posterior draws of \eqn{\sigma} when 
#' the (conjugate) prior on the variance is 
#' \eqn{\sigma^2 ~ Scaled Inv-\chi^2(\alpha, \beta)}, where \eqn{\alpha} 
#' (\code{prior_alpha}) is the shape and \eqn{\beta} (prior_beta) is the scale. 
#' The posterior distribution of the variance parameter is
#' \deqn{
#'   Scaled Inv-\chi^2( m + \alpha, \alpha * \beta * \sum(Y_i - \mu_i)^2 / 
#'     (m + \alpha) ) ,
#' }
#' where \eqn{m} is the number of residuals.
#'
#' @param residuals A vector of residuals \eqn{Y_i - \mu_i}.
#' @param prior_lower,prior_upper Numeric values specifying the lower and upper 
#'   bounds for the uniform prior on \eqn{\sigma}.
#' @param prior_alpha,prior_beta Numeric values representing the shape 
#'   (\code{alpha}) and scale (\code{beta}) of the inverse Chi-squared prior 
#'   distribution on the variance \eqn{\sigma^2}.
#' @param n An integer representing the number of posterior samples to draw 
#'   (defaults to 1).
#'
#' @examples
#' draw_sd_uniform(residuals   = rnorm(130, 0, 2),
#'                 prior_lower = 1,
#'                 prior_upper = 5)
#' 
#' draw_sd_invchisq(residuals   = rnorm(130, 0, 2),
#'                  prior_alpha = 2,
#'                  prior_beta  = 1,
#'                  n = 2)
#'
#' @name draw_sd
NULL


#' @rdname draw_sd
#'
#' @export
draw_sd_uniform <- function(residuals,
                            prior_lower,
                            prior_upper,
                            n = 1) {
  ## Gamma distribution parameters
  shape <- (length(residuals) - 1) / 2
  rate  <- sum(residuals^2) / 2

  ## Bounds for precision (1 / sd^2), given sd \in [prior_lower, prior_upper]
  prec_lower <- 1 / (prior_upper ^ 2)
  prec_upper <- 1 / (prior_lower ^ 2)

  ## Probability that unconstrained precision lies between these bounds
  cdf_upper <- stats::pgamma(prec_upper, shape, rate)
  cdf_lower <- stats::pgamma(prec_lower, shape, rate)
  cdf_delta <- cdf_upper - cdf_lower

  ## Apply inverse distribution method to the constrained distribution
  prec <- stats::qgamma(cdf_lower + stats::runif(n) * cdf_delta, shape, rate)
  1 / sqrt(prec)
}


#' @rdname draw_sd
#' 
#' @export
draw_sd_invchisq <- function(residuals,
                             prior_alpha,
                             prior_beta,
                             n = 1) {
  alpha_posterior <- length(residuals) + prior_alpha
  beta_posterior  <- (prior_alpha * prior_beta + sum(residuals^2)) / 
                        alpha_posterior

  sqrt(rinvchisq(n,
                 alpha_posterior,
                 beta_posterior))
}
