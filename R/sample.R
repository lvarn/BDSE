#' Weighted sample without replacement from table
#'
#' @param n A positive integer; values are sampled from \code{1:n}.
#' @param size Positive integer specifying how many values to sample.
#' @param counts Numeric vector of length \code{n} giving the maximumum number
#'   of each integer that can be drawn.
#' @param probs Numeric vector of length \code{n} giving weights for each
#'   integer.
#'
#' @return A vector of values in \code{1:n} such that \code{i} does not appear
#'   more than \code{count[i]} times.
#'
#' @examples
#'    sample_expj(20,
#'                50,
#'                sample(10, 20, replace = TRUE),
#'                round(runif(20), 2))
#'
#' @export
sample_expj <- function(n,
                        size,
                        counts,
                        probs) {
  if (length(counts) != n) 
    stop("length(counts), ", length(counts), " != n,", n)
  
  tmp <- wrswoR::sample_int_expj(sum(counts),
                                 size = size,
                                 prob = rep(probs, times = counts))
  indices <- rep(seq_len(n), times = counts)
  indices[tmp]
}


#' Draw without replacement
draw <- function(n,
                 size,
                 counts,
                 probs) {
  indices <- sample_expj(n, size, counts, probs)
  make_count(n, indices)
}


#' Weighted sample without replacement with repeated values
#'
#' @param n Number of items to choose from.
#' @param size Number of items to choose.
#' @param counts Vector of length \code{n} giving the number of times each value
#'   occurs.
#' @param probs Vector of weights giving relative likelihood of sampling each
#'   value.
#'
#' @return A vector of length \code{size} with elements of \code{1:n}.
sample_rank <- function (n, 
                         size, 
                         counts, 
                         probs)
{
  if (length(counts) != n)
    counts <- rep(counts, times = rep(ceiling(n / length(counts)),
                                      length(counts)))[1:n]
  indices <- rep(1:n, times = counts)
  indices[head(order(rexp(length(indices)) / rep(probs, times = counts)), size)]
}


#' Sample location
#'
#' @import data.table
sample_location <- function(dt,
                            cluster_col,
                            prob,
                            n_clusters) {
  agree <- stats::rbinom(nrow(dt), 1, prob) == 1
  new_clus <- agree
  new_clus[agree]  <- dt[(agree), eval(as.name(cluster_col))]
  new_clus[!agree] <-
    sample_without_collisions(seq_len(n_clusters),
                              dt[(!agree),
                              eval(as.name(cluster_col))])
  new_clus
}


#' Sample without collision
#'
#' Sample values from support excluding values in bad.
#'
#' @import data.table
sample_without_collisions <- function(support,
                                      bad,
                                      prob = NULL) {
  out <- sample(support, length(bad), replace = TRUE, prob = prob)
  collisions <- out == bad

  while (any(collisions)) {
    out[collisions] <- sample(support, sum(collisions),
                              replace = TRUE, prob = prob)
    collisions <- out == bad
  }
  out
}

