#' Draw from posterior distribution of target population size
#' 
#' Draws a sample from the conditional posterior distribution of the size 
#' (\eqn{N}) of the target population, where
#' \deqn{
#'   \pi(N | data, \psi) \propto N! / (N - k)! 
#'     q(data, \psi)^{N - k}  \pi(N) ,
#' }
#' for different prior distributions \eqn{\pi(N)}, where \eqn{q(data, \psi)} is 
#' the probability of failure (which depends on the data and other parameters), 
#' and \eqn{\psi} is a vector of all other parameters of the model.
#' 
#' \code{draw_n_target_jeffreys} makes draws from the posterior distribution 
#' where the prior distribution on \eqn{N} is the Jeffrey's prior
#' \eqn{p(N) \propto 1 / N}. Here, the posterior distribution of \eqn{N} is a 
#' negative binomial distribution with parameters \code{n_success}, 
#' and \code{p_success} (\eqn{= 1 - q(data, \psi)}).
#' 
#' \code{draw_n_target_uniform} makes draws from the posterior distribution
#' where the prior distribution on \eqn{N} is uniform on (\code{prior_lower},
#' \code{prior_upper}). Here, the posterior distribution is not a standard
#' distribution, and therefore, rejection sampling is used to draw samples from
#' this distribution.
#' 
#' If \code{prior_lower} is set to an integer much larger than \code{n_success},
#' then, when \code{p_success} is large (or equivalently the probability of 
#' failure \eqn{q(data, \psi)} is small), \code{n_try_max} might be reached 
#' before a draw is accepted. 
#' 
#' \code{prior_lower} defaults to \code{n_success}, and \code{prior_upper} 
#' defaults to 
#' \code{n_success/p_success + 3 * (n_success/p_success - n_success)}.
#' 
#' @param n_success An integer representing the number of successes.
#' @param p_success A numeric value in \eqn{(0, 1)} representing the 
#'   probability of success in each Bernoulli trial.
#' @param prior_lower,prior_upper Integer values specifying the minimum 
#'   and maximum of the bounds of the uniform prior distribution on the 
#'   target size. (See Details for defaults.)
#' @param n_try_max An integer specifying the maximum number of tries for the
#'   rejection sampling algorithm (defaults to 10,000).
#' 
#' @return An integer (\eqn{\ge} \code{n_success}).
#' 
#' @examples
#' p_failure <- 0.2
#' 
#' draw_n_target_jeffreys(n_success = 200,
#'                        p_success = 1 - p_failure)
#' 
#' \dontrun{
#' # likely to throw error since p_failure is low and prior_lower is too high
#'   draw_n_target_uniform(n_success   = 200,
#'                         p_success   = 1 - p_failure,
#'                         prior_lower = 500,
#'                         prior_upper = 1000)
#' }
#' 
#' @name draw_n_target
NULL


#' @rdname draw_n_target
#' 
#' @export
draw_n_target_jeffreys <- function(n_success,
                                   p_success,
                                   n_draws = 1) {
  n_failure <- stats::rnbinom(n    = n_draws,
                              size = n_success, 
                              prob = p_success)
  n_success + n_failure
}


#' @rdname draw_n_target
#' 
#' @export
draw_n_target_uniform <- function(n_success,
                                  p_success,
                                  n_try_max   = 10000,
                                  prior_lower = n_success,
                                  prior_upper = n_success/p_success + 
                                    3 * (n_success/p_success - n_success)) {
  if (n_success > prior_lower | prior_lower > prior_upper)
    stop("the order n_success <= prior_lower <= prior_upper must hold")

  n_star            <- n_success / p_success
  log_posterior_max <- calc_loglike_binom_n_trial(N = n_star, 
                                                  k = n_success,
                                                  p = p_success)
 
  for (i in 1:n_try_max) {
    n_target <- prior_lower + sample(prior_upper - prior_lower, size = 1)
    log_posterior_n_target <- calc_loglike_binom_n_trial(N = n_target, 
                                                         k = n_success,
                                                         p = p_success)
    likelihood_ratio <- exp(log_posterior_n_target - log_posterior_max)

    if (stats::runif(1) < likelihood_ratio) break
    if (i == n_try_max)
      stop("maximum number of tries reached without accepting drawn sample")
    }

  n_target
}
