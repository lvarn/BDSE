#' Wrapper around \code{\link{do.call}} allowing unused arguments
#'
#' The \code{args} list can contain values not to be passed to \code{what}.
#'
#' @inheritParams base::do.call
#'
#' @return The result of the (evaluated) function call.
do_call <- function(what,
                    args,
                    quote = FALSE,
                    envir = parent.frame()) {
  args <- args[intersect(names(args),
                         names(formals(what)))]
  do.call(what, args, quote, envir)
}


#' Vectorised is.null
#'
#' @param ... Objects.
#'
#' @return Logical vector.
#'
#' @export
is_null <- function(...) vapply(list(...), is.null, TRUE)


#' Fast tapply
#'
#' @param X A numeric vector.
#' @param INDEX A vector of the same length as \code{X}.
#' @param FUN A function.
#'
#' @return \code{tapply(X, INDEX, FUN = FUN)}
t_apply <- function(X,
                    INDEX,
                    FUN) {
  ord <- order(INDEX)
  fastmatch::ctapply(X[ord], INDEX[ord], FUN = FUN)
}


#' Order elements of a list or vector by name
#'
#' Orders elements of a named list or vector by element name.
#'
#' @param input A named list or vector.
#' @param as_integer Logical, indicating whether the names must be coerced to
#'   integers before ordering.
#'
#' @return The ordered input vector or list.
#'
#' @examples
#' BDSE:::order_by_name(c('2' = 26.4, '1' = 9.5, '3' = 12))
#'
#' BDSE:::order_by_name(list('2' = 26.4, '1' = 9.5, '3' = 12))
order_by_name <- function(input,
                          as_integer = TRUE) {
  if (as_integer)
    input[order(as.integer(names(input)))]
  else
    input[order(names(input))]
}


#' Logarithm of the Normal density function
#' 
#' Calculates the natural logarithm of a vector of normal density values.
#'
#' @param x A vector of numeric values.
#' @param mean A vector of means.
#' @param sd A vector of standard deviations.
#'
#' @return A vector of log densities.
log_dnorm <- function(x,
                      mean,
                      sd) {
  stats::dnorm(x, mean, sd, log = TRUE)
}


#' Draw from a truncated normal
#' 
#' @param n An integer specifying the number of variates to generate.
#' @param mean A vector of means.
#' @param sd A vector of standard deviations.
#' @param min,max The truncation bounds; default to \code{-Inf} and \code{Inf}.
#'
#' @examples 
#' BDSE:::rtruncated_norm(n    = 2, 
#'                        mean = c(0, 10), 
#'                        sd   = c(1, 2.5),
#'                        min  = c(-3, 5), 
#'                        max  = c(3, 25))
#' 
#' @return A vector of length \code{n}.
## Matt's version:
rtruncated_norm <- function(n,
                            mean,
                            sd,
                            min  = -Inf,
                            max  = Inf,
                            ntry = 1e6) {
  out <- stats::rnorm(n, mean, sd)
  i <- 0
  bad_indices <- out < min | out > max
  n_bad <- sum(bad_indices)
  while (n_bad > 0 & i < ntry) {
    out[bad_indices] <- stats::rnorm(n_bad, mean, sd)
    bad_indices <- out < min | out > max
    n_bad <- sum(bad_indices)
    i <- i + 1
  }
  out
}

## TODO: Test this (keep in mind Gibbs L2b)
# rtruncated_norm <- function(n,
#                             mean,
#                             sd,
#                             min = -Inf,
#                             max = Inf) {
#   truncnorm::rtruncnorm(n    = n, 
#                         a    = min,
#                         b    = max,
#                         mean = mean,
#                         sd   = sd)
# }
#

#' Draw from scaled inverse Chi-squared distribution
#'
#' Generates samples from the Scaled Inverse Chi-squared distribution.
#'
#' Samples are drawn from the Gamma distribution with shape-scale
#' parameterization:
#' \deqn{
#'   Gamma( r / 2, 2 / (r * s^2) ) ,
#' }
#' and then transformed. When X has the above Gamma distribution, then
#' \eqn{1 / X} has a scaled inverse Chi-squared distribution with \eqn{r}
#' (\code{deg_free}) degrees of freedom and scale parameter \eqn{s^2}
#' (\code{ssq}).
#'
#' When \code{ssq} is set to \code{1 / deg_free}, then the drawn samples are
#' from an inverse Chi-squared distribution with degrees of freedom
#' \code{deg_free}.
#'
#' @param n An integer specifying the number of observations.
#' @param deg_free An integer specifying the degrees of freedom.
#' @param ssq The scale parameter, \eqn{s^2}.
#'
#' @return A vector of length \code{n}.
rinvchisq <- function(n,
                      deg_free,
                      ssq) {
  1 / stats::rgamma(n, shape = deg_free / 2, rate = (deg_free * ssq) / 2)
}


#' Inverse logit function
#'
#' Given a real value, returns the inverse-logit or logistic function value.
#'
#' The logistic function \eqn{1 / (1 + \exp(-x))} is a mapping from
#' \eqn{(-\Infty, \Infty)} to \eqn{(0, 1)}.
#'
#' If input is not a vector (e.g. is a matrix), function produces a warning,
#' while returning a vector.
#'
#' @param x A vector of numeric values.
#'
#' @return A vector of probabilities.
#'
#' @examples
#' BDSE:::invlogit(matrix(c(1:4), ncol = 2))
#' 
#' ## warns if input is coerced to vector
#' BDSE:::invlogit(matrix(c(1:4), ncol = 1))
#' BDSE:::invlogit(matrix(c(1:4), nrow = 1))
#' 
## TODO: what should the expected behavior be? drop a dimension?
invlogit <- function (x, simplify = TRUE) {
  if (length(dim(x)) == 2L & any(dim(x) == 1L) & simplify) {
    warning("output is coerced to vector; set simplify = FALSE to stop this \n")
    as.vector(1 / (1 + exp(-x)))
  } else
      1 / (1 + exp(-x))
}


#' Calculate medians
#' 
#' Calculates the column/row medians of a matrix or data.frame.
#' 
#' @param x A numeric matrix or data.frame with numeric columns.
#' 
#' @return A vector of medians.
#' 
#' @examples
#' colMedians(matrix(1:20, ncol = 2))
#' colMedians(data.frame(a = 1:10, b = 11:20))
#' 
#' @export
colMedians <- function(x) as.vector(apply(x, 2, median))

#' @rdname colMedians
#' @examples
#' rowMedians(matrix(1:20, ncol = 2))
#' rowMedians(data.frame(a = 1:10, b = 11:20))
#' @export
rowMedians <- function(x) as.vector(apply(x, 1, median))


#' Calculate quantiles
#' 
#' Calculates the column/row quantiles of a matrix or data.frame.
#' 
#' @param x A numeric matrix or data.frame with numeric columns.
#' 
#' @return A vector of quantiles.
#' 
#' @examples
#' colMedians(matrix(1:20, ncol = 2))
#' colQuantiles(matrix(1:20, ncol = 2), p = 0.5)
#' colMedians(data.frame(a = 1:10, b = 11:20))
#' colQuantiles(data.frame(a = 1:10, b = 11:20), p = 0.5)
#' 
#' @export
colQuantiles <- function(x, p) as.vector(apply(x, 2, quantile, prob = p))

#' @rdname colQuantiles
#' @examples
#' rowMedians(matrix(1:20, ncol = 2))
#' rowQuantiles(matrix(1:20, ncol = 2), p = 0.5)
#' rowMedians(data.frame(a = 1:10, b = 11:20))
#' rowQuantiles(data.frame(a = 1:10, b = 11:20), p = 0.5)
#' @export
rowQuantiles <- function(x, p) as.vector(apply(x, 1, quantile, prob = p))


#' Min-max normalize a numeric vector
#'
#' Normalizes a numeric vector to within a specified range (by default
#' \eqn{[-1, 1}).
#'
#' @param var A numeric vector.
#' @param l_min A scalar specifying the lower limit of the normalized range.
#' @param l_max A scalar specifying the upper limit of the normalized range.
#'
#' @return A numeric vector of normalized values.
#'
#' @examples
#' (x <- sample(100, 10))
#' normalize(x)
#' 
#' @export
normalize <- function(var,
                      l_min = -1,
                      l_max = 1) {
  l_min + (l_max - l_min) * (var - min(var)) / (max(var) - min(var))
}


#' Make counts vector
#'
#' @param n Number of elements.
#' @param rows Indices of repeated elements.
#'
#' @return Counts of the number of occurences of each element.
#'
#' @examples
#' BDSE:::make_count(10, c(1, 1, 2, 1, 7))
make_count <- function(n,
                       rows) {
  out <- rep(0, n)
  for (i in rows) {
    out[i] <- out[i] + 1
  }
  out
}

## cpp version
# make_count <- Rcpp::cppFunction(
#   'NumericVector make_count(int n,
#                             NumericVector rows) {
#   NumericVector out(n);
#   for(int i = 0; i < rows.size(); i++) {
#     out[rows[i] - 1]++;
#   }
#   return out;
#   }'
# )


#' Indicator function
#'
#' @param burn_in Positive integer.
#' @param thinning Positive integer. If not provided can be calculated from
#'   \code{count} and \code{nsim}.
#' @param count Positive integer. Only used if \code{thinning} is not provided.
#' @param n_sim Positive integer. Only used if \code{thinning} is not provided.
#'
#' @return A function \eqn{N \rightarrow N} that returns true for every
#'   \code{thinning}th value after the \code{burnin}.
#'
#' @export
record_count <- function(burn_in,
                         thinning,
                         count,
                         n_sim) {
  if (!exists("thinning"))
    thinning <- floor((n_sim - burn_in) / count)
  function(i)
    i > burn_in && (i - burn_in) %% thinning == 0
}


#' Restrict array
#'
#' @param arr Array to restrict in one dimension.
#' @param index Index along dimension \code{dim} to restrict to.
#' @param dim Dimension to collapse.
#'
#' @return An array \code{arr[, ..., index, ...]} where the \code{index} value
#'   occurs in the \code{dim} place.
#'
#' @examples
#'   arr <- array(sample(360, 360),
#'                dim = 3:6)
#'   identical(arr[1, , , ], BDSE:::first(arr))
first <- function(arr,
                  index = 1,
                  dim   = 1) {
  e              <- parent.frame()
  indices        <- rep(list(bquote()), length(dim(arr)))
  indices[[dim]] <- rlang::enexpr(index)
  call           <- as.call(c(list(as.name("["), substitute(arr)),
                              indices))
  eval(call,
       envir = e)
}


#' Assign to first element in first dimension of an array
#'
#' Intended behaviour: \code{arr[1, ...] <- arr1}.
#'
#' \code{arr} and \code{arr1} must exist in the parent environment.
#'
#' @param arr Array to assign values to.
#' @param arr1 A vector, matrix or array containing values to assign to the
#' first element/dimension of \code{arr}.
#'
#' @return Invisibly returns \code{arr1}.
#'
#' @examples
#' # Example 1:
#' a <- array(1, dim = 5)
#' a_1 <- 0
#' BDSE:::assign_first(a, a_1)
#' a
#'
#' # Example 2:
#' a <- array(1:12, dim = c(3, 2, 2))
#' a_1 <- matrix(rep(100, 4), nrow = 2)
#' BDSE:::assign_first(a, a_1)
#' a
assign_first <- function(arr,
                         arr1) {
  arr_call  <- substitute(arr)
  arr1_call <- substitute(arr1)
  e         <- parent.frame()
  call      <- as.call(c(as.name("<-"),
                         as.call(c(
                           list(as.name("["), arr_call),
                           c(1, rep(list(bquote()), length(dim(arr)) - 1)))),
                         arr1_call))
  eval(call, envir = e)
}


#' Zip together lists
#'
#' @param ... Vectors or lists, all of the same length.
#'
#' @return A list with its \eqn{i}th element being the list of the \eqn{i}th
#'   elements of the inputs.
#'
#' @examples
#'   zip(x = 1:3, y = letters[1:3])
#'
#' @export
zip <- function(...) {
  lists <- list(...)
  n     <- unique(vapply(lists, length, 1))
  if (length(n) > 1) stop ("Lists must have the same length")

  if (length(n) == 0)  list()
  else lapply(`names<-`(seq_along(lists[[1]]), names(lists[[1]])),
              function(i) lapply(lists, `[[`, i))
}


#' Make colors for plots
#'
#' @param index_var A vector to color by.
#'
#' @return A vector of colors grouped by \code{index_var}.
#'
#' @export
color_by <- function(index_var,
                     luminance = 70,
                     chroma    = 100) {
  indices   <- unique(index_var)
  hues      <- seq(5, 355, length = length(indices))
  index_var <- match(index_var, indices)
  hcl(h = hues, l = luminance, c = chroma)[index_var]
}


#' Function to make Gibbs sampler result file path
#'
#' This (internal) function makes the file path for the Gibbs samplers.
#'
#' @param path A string specifying the path where the file must be stored.
#' @param model A string specifying the model fit in the sampler; options are
#'   "L1", "L2a" or "L2b". 
#' @param date_time A string specifying the date and time the sample is run.
#' @param seed An integer specifying the random seed used in the sampler.
#' 
#' @return A string, the results file path.
#' 
#' @examples
#' BDSE:::make_file_path("dir/", "L2a", "20120131", 123)
#' BDSE:::make_file_path("dir", "L2a", "20120131", 123)
make_file_path <- function(path, 
                           model, 
                           date_time, 
                           seed) {
  file_name <- paste0(model, "_", date_time, "_seed_", as.character(seed), ".RDS")
  file.path(path, file_name)
}


#' Calculate logistic probabilities
#'
#' Calculates the logistic regression probabilities given a design matrix,
#' coefficients and offsets. This is a wrapper for \code{\link[BDSE]{invlogit}}.
#'
#' @param X A design matrix.
#' @param coefs A vector of coefficients.
#' @param offsets A numeric vector of offsets for the linear component; defaults
#'   to zeros.
#'
#' @return A vector of logistic probabilities with length equal to number of
#'   rows in the design matrix.
#'
#' @examples
#' dataset   <- data.frame(x = sample(10L, 20, replace = TRUE))
#' n_clus    <- 3L
#' cluster   <- sample(n_clus, nrow(dataset), replace = TRUE)
#' clus_effs <- runif(n_clus, -2, 2)
#'
#' calc_probs(X       = model.matrix(~ x, dataset),
#'            coefs   = c(0.5, 0.25),
#'            offsets = clus_effs[cluster])
#'
#' calc_probs(X     = model.matrix(~ x, dataset),
#'            coefs = c(0.5, 0.25))
#'
#' @export
calc_probs <- function(X,
                       coefs,
                       offsets = rep(0, nrow(X))) {
  BDSE:::invlogit(as.vector(X %*% coefs) + offsets)
}


#' Calculate list exclusion probabilities
#'
#' Calculates the probability of being missed by all lists, given the observed
#' covariates data.
#'
#' The probability of not being captured by either of the two lists is \deqn{
#' P{L_1 = 0, L_2 = 0 | X = x, \phi_1, \phi_2} = P{L_1 = 0| X = x, \phi_1} *
#' P{L_2 = 0 | X = x, \phi_2} } assuming that, conditional on the covariate
#' values, inclusion on List 1 and List 2 are independent. For list \eqn{j},
#' \eqn{\phi_j} is the vector of parameters, and includes individual-level
#' effects (\code{coefficients}) and group-level effects (which make up the
#' \code{offsets}).
#'
#' @param ... Named lists with elements: \code{X} (a design matrix of dimension
#'   \eqn{n x p}); \code{coefficients} (a \eqn{p}-dimensional vector of
#'   individual-level effects); and \code{offsets} (a vector of offsets, e.g.
#'   group-level effects for individuals). 
#'
#' @return A vector of probabilities.
#'
#' @examples
#' n          <- 1000  # number of observations
#' n_clusters <- 10    # number of clusters
#' clusters   <- sample(n_clusters, n, replace = TRUE)
#' 
#' X          <- cbind('(Intercept)' = rep(1, n),
#'                     sexM          = rbinom(n, 1, 0.5),
#'                     age           = sample(85, n, replace = TRUE))
#' coeffs1    <- c(0.5, -0.75, 0.01)
#' coeffs2    <- c(0.5, 0.9, 0.05)
#' effects1   <- rnorm(n_clusters, 0, 1)
#' effects2   <- rnorm(n_clusters, 0, 0.5) 
#' offsets1   <- effects1[clusters]
#' offsets2   <- effects2[clusters]
#' 
#' ## Example 1:
#' p_mis <- BDSE:::calc_missing_probabilities(
#'            list(X            = X,
#'                 coefficients = coeffs1,
#'                 offsets      = offsets1),
#'            list(X            = X,
#'                 coefficients = coeffs2,
#'                 offsets      = offsets2))
#' head(p_mis)
#' 
#' ## Example 2: will have two unique probabilities
#' unique(BDSE:::calc_missing_probabilities(
#'          list(X            = X,
#'               coefficients = c(0.5, 1, 0)),
#'          list(X            = X,
#'               coefficients = c(0.1, 0.5, 0))))
calc_missing_probabilities <- function(...) {
  exclusion_probabilities <-
    lapply(list(...),
           function(x) {
             offsets <- if (is.null(x[["offsets"]])) 0 else x[["offsets"]]
             1 - calc_probs(x[["X"]], x[["coefficients"]], offsets)
           })
  
  Reduce(`*`, exclusion_probabilities)
}


