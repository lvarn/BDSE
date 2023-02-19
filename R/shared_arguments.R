#' Shared arguments
#'
#' @param X A numeric (design) matrix with dimension \eqn{n} x \eqn{p}, where
#'   \code{n} is the number of individuals, and \eqn{p} is the number of
#'   covariates (including the intercept).
#' @param y A binary (0, 1) vector of length \eqn{n}, representing observations
#'   of the response variable.
#' @param prob A vector of (success) probabilities having the same length as the
#'   reposnse vector.
#' @param p_cov The vector of coverage probabilities; defaults to \code{NULL}.
#' @param p_sens The vector of linkage sensitivities; defaults to \code{NULL}.
#' @param cluster_effects1,cluster_effects2 Optional argument giving the cluster
#'   effects used when calculating coverage probabilities.
#' @param coefficients A numeric vector of length \eqn{p + 1} of
#'   individual-level effects in the logistic regression model, the first
#'   element of which is the intercept.
#' @param cluster_effects A numeric vector of length \code{n_clusters} of
#'   group-level (cluster) effects, where \code{n_clusters} is the number of
#'   unique clusters in the population.
#' @param clusters A vector of length \eqn{n}, with values in
#'   \code{1:n_clusters}, representing the cluster membership of the individual
#'   observations.
#' @param n_clusters An integer specifying the number of unique clusters in the
#'   population, some of which may not be represented in the observed data.
#' @param offsets An optional numeric vector of length \eqn{n} specifying the
#'   values added to the linear component to account for effects not dependent
#'   on \code{X}.
#' @param weights An optional numeric vector of length \eqn{n} specifying the
#'   sample weights.
#' @param counts A numeric vector of counts of the number of observations in
#'   each corresponding row of the design matrix \code{X}.
#' @param formula1,formula2 Each a formula for creating the design matrix for
#'   the coverage model of the corresponding list (list 1 or 2).
#' @param ntry An integer specifying the maximum number of tries until \code{n}
#'   truncated normal variates are generated.
#' @param periodic_save A named list with two elemnts: \code{every}, an integer
#'   \eqn{n} specifying that every \eqn{n}th iteration must be saved, and
#'   \code{where}, a string specifying the directory where the results should be
#'   saved. (Defaults to \code{NULL}.)
#' @param print_every An integer \eqn{n} specifying that after every \eqn{n}th
#'   iteration a progress report should be printed (to console).
#' @param record_count A function indicating which counts (iterations) should be
#'   stored in the output object; see \code{\link[BDSE]{record_count}}.
#' @param return_obs_dist Logical, indicating whether the observed frequency
#'   distribution should be stored in the output object; defaults to
#'   \code{TRUE}.
#'
#' @name arguments
NULL
