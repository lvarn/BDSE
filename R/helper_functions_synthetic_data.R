#' Get names of dummified data.frame
#'
#' This function returns the names of the input data.frame after 'dummifying'
#' categorical variables.
#'
#' @param dataset A data.frame.
#' @param model_formula A formula, which will be passed to \code{model.matrix()}
#'   in the function to dummify the input data.frame. (Defaults to \code{~ .} if
#'   not specified.)
#'
#' @return A single-column character matrix of dummified variable names.
#'
#' @examples
#' n <- 1000
#' X_df <- data.frame(income    = runif(n, 1, 1000),
#'                    sex       = sample(c("F", "M"), n, replace = TRUE),
#'                    age       = sample(100, n, replace = TRUE),
#'                    ethnicity = sample(LETTERS[1:5], n, replace = TRUE))
#'
#' get_dummified_names(dataset = X_df)
#'
#' get_dummified_names(dataset       = X_df,
#'                     model_formula = ~ age + ethnicity)
#'
#' @export
get_dummified_names <- function(dataset,
                                model_formula = ~ .) {
   as.matrix(
     colnames(
       stats::model.matrix(model_formula, dataset)))
}


#' Extend a vector of coefficients
#'
#' This function, given a design matrix \code{X} and a vector of coefficients
#' \code{b}, extends \code{b} with 0-valued elements to be compatible with
#' \code{X}, so that \code{X \%*\% b} is correctly computed.
#'
#' This function is a helper function for \code{generate_target_dataset}. It is
#' used to extend the coefficients vector to match the number of columns in the
#' design matrix, when only some coefficients are specified (see examples).
#'
#' @param X A design matrix with column names.
#' @param b A named numeric vector, with names matching column names of
#'   \code{X}.
#'
#' @return A named numeric vector of length equal to \code{ncol(X)}.
#'
#' @examples
#' n       <- 50
#' dataset <- data.frame(age_grp   = factor(sample(seq(10, 65, 5), n,
#'                                                 replace = TRUE)),
#'                       income    = runif(n, 1, 50000),
#'                       sex       = sample(c("M", "F"), n, replace = TRUE),
#'                       hours_grp = sample(c("[0,10)", "[10,40)", "[40,Inf)"),
#'                                          n, replace = TRUE))
#' X <- model.matrix(~ ., data = dataset)
#' head(X)
#'
#' b <- c('(Intercept)'      = 0.5,
#'        age_grp15          = 2,
#'        age_grp35          = 1.5,
#'        'hours_grp[10,40)' = 3)
#'
#' b_extended <- extend_coefficients(X = X, b = b)
#'
#' as.matrix(b)
#' as.matrix(b_extended)
#' 
#' @export
extend_coefficients <- function(X, b) {
  ## Check that names(b) is a strict subset of colnames(X)
  names_X <- colnames(X)
  names_b <- names(b)

  if(!all(names_b %in% names_X))
    stop("coefficient (vector) names must be a subset of the column names of ",
         "the corresponding design matrix")

  if("(Intercept)" %in% names_X & !"(Intercept)" %in% names_b)
    warning("the coefficients vector does not include an intercept term: ",
            "the value for '(Intercept)' will be set to 0")

  ## Complete the coefficients vector by replacing unspecified terms with 0s
  b_extended <- `names<-`(rep(0, length(names_X)), names_X)
  
  b_extended[match(names_b, names(b_extended))] <- b

  b_extended
}

