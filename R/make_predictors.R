## TODO: write generic functions for making ethnicity interactions

#' Make ethnicity indicators
#'
#' This function converts the vector of ethnicity codes into a data.frame of
#' ethnicity indicators for single ethnicities.
#'
#' In \code{base_dataset}, the \code{ethnicity} column is the integer-valued
#' ethnicity codes, where single digits refer to single-ethnicity groups and 
#' double digits refer to bi-ethnicity groups. This function converts this
#' ethnicity encoding to a data.frame of single-ethnicity indicators, where each
#' individual (row) has at least one single-ethnicity indicator set to 1.
#'
#' @param ethnicity A vector of numeric ethnicity codes.
#' @param single_ethnicities A vector specifying the values of the
#'   unique single-ethnicity levels. (Defaults to the single ethnicity codes
#'   \code{c(1:6, 9)}.)
#'
#' @return A data.frame of ethnicity indicators.
#'
#' @examples
#' dataset <- BDSE::base_dataset
#' table(dataset$ethnicity)
#' 
#' ethnicity_df <- make_ethnicity_indicators(dataset[["ethnicity"]])
#' 
#' head(ethnicity_df)
#' colSums(ethnicity_df)
#' 
#' ## Checks:
#' all(rowSums(ethnicity_df) > 0)
#' all(ethnicity_df[ethnicity_df[, 7] == 1, -7] == 0)  # 7 = ethnicity9 = "other"
#' 
#' @export
make_ethnicity_indicators <- function(ethnicity,
                                      single_ethnicities = c(1:6, 9)) {
  cols <- lapply(`names<-`(single_ethnicities, 
                           paste0("ethnicity", single_ethnicities)),
                 function(i) as.integer(grepl(i, ethnicity)))
  
  as.data.frame(cols)
}


#' Make ethnicity interactions
#' 
#' This functions creates interaction indicators from the single-ethnicity 
#' indicators (which are assumed to exist in the dataset).
#' 
#' @param interactions A list, each element of which is a 2-element list with
#'   its first element the code for the single indicator, and its second
#'   element, the codes for the single ethnicities it interacts with. (See
#'   defaults.)
#' @param base_name A string represnting the base name of the variable (defaults
#'   to \code{"ethnicity"}).
#' @param dataset A data.frame with the single ethnicity indicators as columns.
#' 
#' @return A data.frame of ethnicity interactions indicators.
#' 
#' @examples 
#' dataset <- BDSE::base_dataset
#' dataset <- cbind(dataset, 
#'                  make_ethnicity_indicators(dataset[["ethnicity"]]))
#' 
#' ethnicity_interactions_df <- make_ethnicity_interactions(dataset = dataset)
#' 
#' head(ethnicity_interactions_df)
#' 
#' @export
make_ethnicity_interactions <- function(dataset,
                                        interactions = 
                                          list(list(1, 2:6),
                                               list(2, 3:6),
                                               list(3, 4:6),
                                               list(4, 5:6),
                                               list(5, 6)), 
                                        base_name    = "ethnicity") {
  matrices <- lapply(interactions,
                     function(x) {
                       M <- dataset[[paste0(base_name, x[[1]])]] * 
                              dataset[paste0(base_name, x[[2]])]
                       `colnames<-`(M, paste0(base_name, x[[1]], ".", x[[2]]))
                     })
  
  do.call(cbind, matrices)
}


#' Create a basis matrix for a polynomial spline
#' 
#' This function is a wrapper for \code{\link[splines]{bs}}, and is used to
#' allow specifying the piecewise polynomial terms of the associated predictor
#' variable in the logistic regression coverage model.
#' 
#' @param variable A numeric vector: the predictor variable the basis matrix is
#'   created from.
#' @param base_name A string specifying the name of the variable to use as
#'   the base name for the columns of the basis matrix. (For example, for "age"
#'   and 10 terms, the column names will be "age1", "age2", ..., "age10".)
#' @param bs_params A named list of parameters passed to
#'   \code{\link[splines]{bs}}.
#' 
#' @return A matrix.
#' 
#' @examples
#' dataset <- data.frame(sex = sample(2, 1000, replace = TRUE),
#'                       age = sample(100, 1000, replace = TRUE))
#' 
#' bs_matrix <- make_basis_matrix(variable    = "age", 
#'                                dataset     = dataset,
#'                                base_name   = "age_bs",
#'                                smooth_type = "bs",
#'                                params      = list(knots = seq(15, 85, 10)))
#' 
#' head(bs_matrix, 10)
#' 
#' @export
# TODO: tidy up these functions (especially .ti).
make_basis_matrix <- function(variables,
                              dataset,
                              base_name,
                              smooth_type,
                              params        = NULL,
                              include_mains = TRUE) {
  switch(smooth_type,
         bs = do.call(make_basis_matrix.bs, 
                      args = list(variable  = unlist(variables),
                                  dataset   = dataset,
                                  base_name = base_name,
                                  bs_params = params)),
         ti = do.call(make_basis_matrix.ti,
                      args = list(variables     = variables,
                                  dataset       = dataset,
                                  base_name     = base_name,
                                  ti_params     = params,
                                  include_mains = include_mains)))
}


#' @rdname make_basis_matrix
#' 
#' @export
make_basis_matrix.bs <- function(variable,
                                 dataset,
                                 base_name,
                                 bs_params) {
  basis_matrix <- do.call(splines::bs, 
                          args = c(list(x = dataset[[variable]]), bs_params))
  
  `colnames<-`(basis_matrix,
               paste0(base_name, 1:ncol(basis_matrix)))
}


#' @rdname make_basis_matrix
#' 
#' @export
make_basis_matrix.ti <- function(variables,
                                 dataset,
                                 base_name,
                                 ti_params,
                                 include_mains = FALSE) {
 
  S <- eval(as.call(c(mgcv::smoothCon,
                    as.call(c(mgcv::ti, lapply(variables, as.name), ti_params)), 
                    data = bquote(dataset))))
 
  basis_matrix <- S[[1]][["X"]] 
  colnames(basis_matrix) <- paste0(base_name, "_i", 1:ncol(basis_matrix))
  
  if (include_mains) {
      m         <- length(S[[1]][["margin"]])
      marginals <- lapply(1:m, 
                          function(i) {
                            X <- S[[1]][["margin"]][[i]][["X"]]
                            `colnames<-`(X, paste0(base_name, "_", 
                                                   i, ".", 1:ncol(X)))
                          })
  
      basis_matrix <- do.call(cbind, c(marginals, list(basis_matrix)))
  }
  
  basis_matrix
}

