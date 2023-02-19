#' A wrapper to generate all errors
#' 
#' TODO: Fill in parameter details and description later
#' 
#' @export
make_synthetic_data <- function(seed,
                                dataset, # target or union
                                coefficients1,
                                coefficients2,
                                formula1,
                                formula2,
                                clusters_sd1,
                                clusters_sd2,
                                clustering_var,
                                formula_sens          = NULL,
                                formula_target        = NULL,
                                clusters_sd_target    = NULL,
                                coefficients_sens     = NULL,
                                coefficients_target   = NULL,
                                clustering_var_target = NULL,
                                n_validated           = 50000 # not used if sens NULL
                                # 
                                ) {
  ## Apply overcoverage to data
  if (!is.null(formula_target))
    dataset <- add_overcoverage(dataset             = dataset,
                                formula_target      = formula_target,
                                coefficients_target = coefficients_target,
                                clustering_var      = clustering_var_target,
                                clusters_sd_target  = clusters_sd_target,
                                seed                = seed)
 
  ## Apply undercoverage to data
  dataset <- add_undercoverage(dataset          = dataset,
                               formula1         = formula1,
                               formula2         = formula2,
                               coefficients1    = coefficients1,
                               coefficients2    = coefficients2,
                               clustering_var   = clustering_var,
                               clusters_sd1     = clusters_sd1,
                               clusters_sd2     = clusters_sd2,
                               seed             = seed,
                               has_overcoverage = !is.null(formula_target))
  
  ## Add linkage error to observed union data
  if (!is.null(formula_sens)) {
    dataset_obs    <- dataset[!dataset[["Y00"]], ]
    dataset_obs_LE <- add_linkage_error(dataset           = dataset_obs,
                                        formula_sens      = formula_sens,
                                        coefficients_sens = coefficients_sens,
                                        seed              = seed,
                                        n_validated       = n_validated)
  }

  list(D        = dataset,
       D_obs_LE = if(!is.null(formula_sens)) dataset_obs_LE else NULL)
}




