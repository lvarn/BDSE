## TODO: get rid of ones we don't need (coda/mcmc already have these functions)
## Add examples if keeping

#' Function to import Gibbs sampler result files
#'
#' Given the directory and file names (corresponding to each chain), reads in	
#' the results files and merges them into a single list.
#'
#' @param directory A string specifying the directory in which the result files
#'   are stored.
#' @param file_names A character vector of file names (with extensions, e.g.
#'   "gibbs_chain1.RDS")
#'
#' @return A list of length equal to the number of chains (files).
#'   
#' @export
make_results_list <- function(directory,
                              file_names) {
  L <- lapply(file_names, 
              function(file)
                readRDS(file.path(directory, file)))
  `names<-`(L, paste0("chain", seq_along(L)))
}


#' Makes summary of Gibbs sampler stored object
#' 
#' Makes a summary data.frame of a Gibbs sampler output (single chain), i.e the	
#' stored object returned by the sampler function \code{run_Gibbs_sampler}.
#'
#' @param result A list with named elements.
#'
#' @return A summary data.frame of names, dimensions, and classes of the input
#'   list.
#'   
#' @export
summarize_single_chain <- function(result) { 
  dims <- lapply(result,
                 function(x) {
                   if (is.array(x) | is.data.frame(x))
                     dim(x)
                   else
                     length(x)
                 })

  property <- 
    data.frame(name  = names(result),
               dim   = vapply(dims, 
                              function(e) paste0(e, collapse = ","),
                              "a"),
               n_dim = vapply(dims, length, 1L),
               list  = vapply(names(result), 
                              function(e) {
                                if (grepl("1$", e)) 
                                  1L
                                else if (grepl("2$", e))
                                  2L
                                else NA_integer_
                              }, 
                              1L),
               class = vapply(result[names(result)], 
                              function(x) class(x)[1], 
                              "a"),
               size  = vapply(result[names(result)], object.size, 1))
  
  `rownames<-`(property, NULL)
}


#' Binds elements of lists
#'
#' Given a list of sampler chain outputs, this function binds the elements of	
#' each list (i.e. chain) into a single matrix, and returns a list of matrices,	
#' where each matrix corresponds to a chain output.
#'
#' @param results_list A list of lists (i.e. the list of chains from the Gibbs
#'   sampler).
#' @param pars A character vector specifying which groups of parameters must be
#'   retained (e.g. "n_target", "coefficients1", "cluster_sds1").
#'
#' @return A list of matrices, where each element of the list coresponds to a 
#' single chain.
#'
#' @export
bind_chains <- function(results_list, 
                        pars) {
  res <- lapply(results_list, 
           function(chain) {
             chain <- chain[grepl(paste0(pars, collapse = "|"), names(chain))]
             chain <- lapply(pars,
                        function(par) {
                          col_names <- 
                            if (length(dim(chain[[par]])) < 2) 
                              par else {
                                if ( is.null(colnames(chain[[par]])) )
                                  paste(par, seq_len(ncol(chain[[par]])), sep = ".") else
                                  paste(par, colnames(chain[[par]]), sep = ".")
                              }
                          `colnames<-`(as.matrix(chain[[par]]), col_names)
                        })
             
             do.call(cbind, chain)
         })
  
  res
}


#' Convert list of matrices to single data.frame
#' 
#' @param results_list A list of matrices (sampler outputs), where each element
#'   of the list corresponds to a chain; a result of a call to
#'   \code{\link[BDSE]{bind_chains}}.
#' @param n_burnin An integer specifying the number of initial iterations to
#'   discard from each chain.
#' @param n_thin An integer specifying the frequency of post-burnin iterations
#'   to keep from each chain.
#'
#' @return A data.frame of thinned post-burn-in posterior samples. 
#' 
#' @export
convert_to_df <- function(results_list, 
                          n_burnin = 0,
                          n_thin   = 1) {
  results_list <- lapply(seq_along(results_list),
                         function(i) {
                           df   <- results_list[[i]]
                           ids  <- seq(n_burnin + 1, nrow(df), n_thin)
                           df   <- df[ids, ]
                           init <- data.frame(chain = i, 
                                              iteration = 1:nrow(df))
                           cbind(init, df) 
                        })

  do.call(rbind, results_list)
}



#' Convert list of matrices to single 3D array
#' 
#' @inheritParams convert_to_df
#' 
#' @return A 3-dimensional array of thinned post-burn-in posterios samples, with
#'   dimensions iterations, chains, and sample.
#' 
#' @export
convert_to_array <- function(results_list, 
                             n_burnin = 0,
                             n_thin   = 1) {
  results_list <- lapply(results_list,
                         function(mat) {
                           ids  <- seq(n_burnin + 1, nrow(mat), n_thin)
                           mat  <- mat[ids, ]
                           # init <- data.frame(iteration = seq_len(nrow(mat)))
                           # mat  <- cbind(as.matrix(init), mat) 
                        })
  
  results_array <- simplify2array(results_list)
  
  aperm(results_array, c(1, 3, 2))
}


#' Get positions of coefficients
#' 
#' Gets the positions of the coefficients from the Gibbs sampelr results for the
#' given list.
#'
#' @param list_id An integer, the list ID (1 or 2).
#' @param true_names A character vector of names of the coefficients appearing
#'   in the order needed.
#' @param pars A character vector of Gibbs sample coefficient/parameter names to
#'   search in.
#' 
#' @return An integer vector of indices.
#' 
get_coefs_idx <- function(list_id,
                          true_names,
                          pars) {
  pattern <- paste(paste0("coefficients", list_id, "\\."), 
                   paste0("coefficients_clus", list_id, "\\."),
                   sep = "|")
  
  idx <- grep(pattern, pars)
  
  gibbs_names <- gsub(pattern, "", pars[idx])
  
  idx[match(true_names, gibbs_names)]
}


#' Extract Gibbs posterior sample for coefficients
#'
#' Extracts the Gibbs sampler samples for the coefficients of the given list ID.
#'
#' @inheritParams get_coefs_idx
#' @param results A data.frame of Gibbs sampler results (the output of
#'   \code{\link[BDSE]{convert_to_df}}).
#' 
#' @return A data.frame of posterior samples.
#' 
#' @export
extract_coefs <- function(list_id,
                          results, 
                          true_names) {
  idx <- BDSE:::get_coefs_idx(list_id = list_id,
                              true_names = true_names,
                              pars = colnames(results))
  
  `colnames<-`(results[, idx], true_names)
}


#' Extract Gibbs posterior sample for cluster effects
#'
#' Extracts the Gibbs sampler samples for the cluster/area effects of the given
#' list ID.
#'
#' @inheritParams get_coefs_idx
#' @param results A data.frame of Gibbs sampler results (the output of
#'   \code{\link[BDSE]{convert_to_df}}).
#' @param clus_names An optional character vector of names for the cluster
#'   effects.
#'
#' @return A data.frame of posterior samples.
#'
#' @export
extract_clus_effects <- function(list_id,
                                 results,
                                 clus_names = NULL) {
  pattern <- paste0("cluster_effects", list_id)
  idx     <- grep(pattern, colnames(results))

  clus_names <- if ( !is.null(clus_names) ) clus_names else 
                  as.character(seq_along(idx))
  
  `colnames<-`(results[, idx], clus_names)
}



#' Extract posterior counts
#' 
#' Extracts the posterior counts vectors for the covariate combinations.
#'  
#' @param results_list A named list of lists, each elemnt of which corresponds
#'   to a single chain of the Gibbs sampler.
#' @param rename_cols A logical scalar indicating whether the count columns
#'   should be renamed (\code{1:n}); defaults to \code{FALSE}. (When
#'   \code{FALSE}, the column names are of the form \code{chainx.counts_xxx},
#'   where \code{x} is a number.)
#' 
#' @return A data.frame of posterior count columns.
#'   
#' @export
get_counts <- function(results_list,
                       rename_cols = FALSE,
                       n_burnin = 0,
                       n_thin   = 1,
                       as_array = FALSE) {
  results <- lapply(results_list, 
                    function(chain) {
                      df <- as.data.frame(chain[["out_counts"]][,
                        grepl("^count", names(chain[["out_counts"]]))
                      ])
                      df[, seq(n_burnin + 1, ncol(df), n_thin)]
                    })
  ## Array or data.frame?
  if (as_array) {
    results <- aperm(simplify2array(lapply(results, t)), c(1, 3, 2))
    dimnames(results)[[3]] <- paste0("theta", seq_len(dim(results)[3]))
    
    if (rename_cols)
      dimnames(results)[[1]] <- seq_len(dim(results)[1])
  } else {
      results <- do.call(cbind, results) # as.data.frame(t(do.call(cbind, results)))

      if (rename_cols)
        colnames(results) <- paste0("counts_", seq_len(ncol(results)))
  }
  
  results
}


#' Calculate credible intervals
#'
#' Calculates either the equal-tail (ETI) or highest-density (HDI) posterior
#' intervals using the credible interval (CI) function
#' \code{\link[bayestestR]{ci}}. An option is provided here to manually set
#' probabilities and get the desired quantiles for the CI. (This option is
#' equivalent to choosing the CI method "ETI", but median values might vary
#' slightly due to implementation here.)
#'
#' If \code{ci_probs} is specified, it is used to get the quantiles for the CI;
#' otherwise, \code{ci_coverage} and \code{ci_method} will be used.
#'
#' @param posterior_sample A data.frame of posterior samples.
#' @param true_values An optional numeric vector specifying true parameter
#'   values.
#' @param mle_value An optional numeric vector specifying the maximum likelihood
#'   estimates.
#' @param ci_coverage A scalar in the interval \eqn{[0, 1]} specifying the
#'   coverage of the CI.
#' @param ci_method A string (either "ETI" or "HDI") specifying the CI method to
#'   be used.
#' @param ci_probs A two-element numeric vector specifying the lower and upper
#'   probabilities for the credible intervals. (Defaults to \code{NULL}).
#'
#' @return A data.frame, where rows correspond to parameters, and the columns
#'   are the lower and higher bound of the CI, the posterior median and the true
#'   values and ML estimates (if provided as input).
#' 
#' @export
## TODO: test and complete this!
calc_ci <- function(posterior_sample,
                    true_values = NULL,
                    mle_values  = NULL, 
                    ci_coverage = 0.90,
                    ci_method   = c("ETI", "HDI"),
                    ci_probs    = NULL) {
  ## ETI is same as quantile(), but will vary slightly due to implementation here
  if ( !is.null(ci_probs) ) {
    ci_mat <- apply(posterior_sample, 2, 
                    function(col) quantile(col, probs = ci_probs))
    ci_df <- `names<-`(as.data.frame(t(ci_mat)), c("CI_low", "CI_high"))
    medians <- lapply(posterior_sample, median)
  } else {
    ## Primarily for HDI
    ci_list <- apply(posterior_sample, 2, 
                     function(x) ci(x, method = ci_method, ci = ci_coverage))
    
    medians <- lapply(names(ci_list),
                      function(par) {
                        vals <- posterior_sample[[par]]
                        idx  <- vals >= ci_list[[par]][["CI_low"]] &
                                vals <= ci_list[[par]][["CI_high"]]
                        
                        median(vals[idx])
                      })
    
    ci_df <- as.data.frame(do.call(rbind, ci_list))[, -1]
  }

  cbind(ci_df, 
        median = do.call(rbind, medians),
        mle    = mle_values,
        true   = true_values)
}


