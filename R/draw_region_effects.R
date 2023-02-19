#' Draw region effects from its conditional posterior distribution
#'
#' This function draws a vector of high-level region effects from its
#' conditional posterior distribution, given a vector of low-level cluster
#' effects. (In this setting, each region contains one or more clusters).
#' 
#' This function generates a draw from the conditional posterior of the means of
#' the normal data model with known standard deviation:
#' \deqn{
#'    Y_j ~ Normal(\mu_{reg(j)} + \eta_j, \sigma^2),
#'  }
#' where \eqn{Y_j} is the \eqn{j}th cluster effect, \eqn{\mu_{reg(j)}} is the
#' effect of the region to which cluster \eqn{j} belongs, and \eqn{\eta_j} is
#' some known offset for cluster \eqn{j} (e.g. a linear combination of
#' cluster-level covariates).
#' 
#' The conjugate prior distribution of the region effects 
#' (\eqn{\mu_1, \dots, \mu_{n_regions}}) is 
#' \preformatted{
#'   Normal(prior_mean, prior_sd).
#' }
#'  
#' @param cluster_effects A numeric vector of low-level cluster effects.
#' @param regions An integer or numeric vector, the same length as
#'   \code{cluster_effects}, specifying the indices of the regions to which the
#'   clusters corresponding to \code{cluster_effects} belong.
#' @param n_regions An integer specifying the number of regions in the
#'   population.
#' @param clusters_sd A numeric scalar representing the between-cluster standard
#'   deviation.
#' @param prior_mean A numeric vector of length \code{n_regions} specifying the
#'   mean of the prior distribution of the region-level effects.
#' @param prior_sd A numeric scalar specifying the standard deviation of the
#'   prior distribution of the region-level effects.
#' @param offsets An optional numeric vector, having the same length as
#'   \code{cluster_effects}, of known offsets (defaults to 0); see Details.
#'
#' @return A numeric vector of length \code{n_regions} of the drawn sample of 
#'   region effects.
#' 
#' @examples 
#' n_clusters      <- 1000     # number of clusters
#' n_regions       <- 20       # number of regions
#' 
#' # draw region indices (allow unobserved regions) 
#' regions         <- sample(sample(n_regions, floor(0.8 * n_regions)),
#'                           size = n_clusters,
#'                           replace = TRUE)
#' table(regions)
#' 
#' cluster_effects <- rnorm(n_clusters, 2, 1)
#' prior_mean      <- rnorm(n_regions)
#' prior_sd        <- 3
#' clusters_sd     <- 1.2
#' 
#' draw_region_effects(cluster_effects = cluster_effects,
#'                     regions         = regions,
#'                     n_regions       = n_regions,
#'                     clusters_sd     = clusters_sd,
#'                     prior_mean      = prior_mean,
#'                     prior_sd        = prior_sd)
#' 
#' @export 
draw_region_effects <-  function(cluster_effects,
                                 regions,       
                                 n_regions,
                                 clusters_sd, 
                                 prior_mean, 
                                 prior_sd,
                                 offsets = 0) {
  if(1L < length(offsets) & length(offsets) < length(cluster_effects))
    stop("'offsets' must be a scalar or a vector of length ",
         "length(cluster_effects)")

  resids <- cluster_effects - offsets  # intercepts only (drop covariates part)
  
  ## Group the cluster effects by region and compute sums and lengths
  resids_sums    <- t_apply(resids, regions, sum)
  resids_lengths <- t_apply(resids, regions, length)

  ## Some regions may not be observed; assign 0 to their mean and length
  n_obs_regions <- length(unique(regions))
  if (n_obs_regions < n_regions) {
      resids_sums_mis        <- rep(0, n_regions - n_obs_regions)
      names(resids_sums_mis) <- setdiff(1:n_regions, unique(regions))
      resids_sums            <- order_by_name(c(resids_sums, 
                                                resids_sums_mis))
      
      resids_lengths_mis        <- rep(0L, n_regions - n_obs_regions)
      names(resids_lengths_mis) <- setdiff(1:n_regions, unique(regions))
      resids_lengths            <- order_by_name(c(resids_lengths, 
                                                   resids_lengths_mis))
  }
  
  ## Compute posterior varaince (1 / precision) and mean for each region
  posterior_precs <- resids_lengths / clusters_sd^2 + 
                       rep(1 / prior_sd^2, n_regions)
  posterior_vars  <- 1 / posterior_precs
  posterior_means <- posterior_vars * (resids_sums / clusters_sd^2 + 
                                         prior_mean / prior_sd^2) 
  
  list(region_effects  = stats::rnorm(n_regions, 
                                      posterior_means, 
                                      sqrt(posterior_vars)),
       posterior_means = posterior_means,
       posterior_sds   = sqrt(posterior_vars))
}

