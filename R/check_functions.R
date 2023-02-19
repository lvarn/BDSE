#' Check overcoverage input
#'
#' Internal function to check for consistency of (List 2) overcoverage
#' information.
#'
#' @inheritParams add_undercoverage
#'
#' @examples
#' BDSE:::check_in_target(FALSE, data.frame(X = 10))
#' BDSE:::check_in_target(TRUE,  data.frame(X = 10))
#' BDSE:::check_in_target(FALSE, data.frame(X = 10, in_target = T))
#' BDSE:::check_in_target(TRUE,  data.frame(X = 10, in_target = T))
#'  
check_in_target <- function(has_overcoverage, 
                            dataset) {
  if ( !has_overcoverage & ("in_target" %in% names(dataset)) )
      warning("'dataset' might include list (L2) overcoverage as indicated by 
              'in_target' column")
  
  if ( has_overcoverage & !("in_target" %in% names(dataset)) )
      stop("'dataset' must include the target inclusion indicator 'in_target' 
           in the presence of list overcoverage")
}


#' Check partial coverage input
#'
#' Internal function to check for consistency of (List 1) partial coverage
#' information.
#'
#' @inheritParams add_undercoverage
#'
#' @examples
#' BDSE:::check_in_list(FALSE, data.frame(X = 10))
#' BDSE:::check_in_list(TRUE,  data.frame(X = 10))
#' BDSE:::check_in_list(FALSE, data.frame(X = 10, in_list = T))
#' BDSE:::check_in_list(TRUE,  data.frame(X = 10, in_list = T))
#'  
check_in_list <- function(has_partial_coverage, 
                          dataset) {
  if ( !has_partial_coverage & ("in_list" %in% names(dataset)) )
      warning("'dataset' might include list (L1) partial coverage by design as 
              indicated by 'in_list' column")
  
  if ( has_partial_coverage & !("in_list" %in% names(dataset)) )
      stop("'dataset' must include the coverage indicator 'in_list' in the
           presence of partial coverage by design")
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## TODO: get rid of ones I don't need!
## this script includes "check" for function arguments;


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to check that a list argument has valid element names
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_list <- function(input_list, 
                       input_list_name) {
  ## 'input_list' must be a list
  if(any(!is.list(input_list), is.data.frame(input_list)))
    stop(gettextf("'%s' must be a list", input_list_name))
  
  ## it cannot have length 0
  if(identical(length(input_list), 0L))
    stop(gettextf("'%s' must have length greater than 0L", input_list_name))

  ## it must have valid element names
  list_names <- names(input_list)
  if(is.null(list_names))
    stop(gettextf("names(%s) cannot be NULL", input_list_name))
  
  if(any(is.na(list_names), !nzchar(list_names)))
    stop(gettextf("names(%s) cannot have blank or missing values", 
                  input_list_name))
  
  if(any(duplicated(list_names)))
    stop(gettextf("names(%s) cannot have duplicates", input_list_name))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to check that a data.frame argument has valid column names
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_dataframe <- function(dataset, dataset_name) {
  if(!is.data.frame(dataset))
    stop(gettextf("'%s' must be a data.frame", dataset_name))
    
  col_names <- names(dataset)
  if(is.null(col_names))
    stop(gettextf("names(%s) cannot be NULL", dataset_name))
  
  if(any(is.na(col_names), !nzchar(col_names)))
    stop(gettextf("names(%s) cannot have blank or missing values", 
                  dataset_name))
  
  if(any(duplicated(col_names)))
      stop(gettextf("names(%s) cannot have duplicates", dataset_name))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to check that a data.frame includes the necessary set of variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_variables_in_dataframe <- function(variable_names,
                                         dataset, dataset_name) {
  check_dataframe(dataset = dataset, dataset_name = dataset_name)
  
  dataset_variables <- names(dataset)  ## variable names in dataset
  for(name in variable_names) {
    if(!name %in% dataset_variables)
      stop(gettextf("'%s' does not include variable called \"%s\"",
                    dataset_name, name))
  }
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to check a numeric variable in a given data.frame
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_numeric_variable <- function(variable_name,
                                   dataset, dataset_name) {
  if(!variable_name %in% names(dataset))
    stop(gettextf("dataset '%s' does not have variable called '%s'",
                  dataset_name, variable_name))

  if(!is.numeric(dataset[[variable_name]]))
    stop(gettextf("'%s' must be of type \"integer\" or \"double\"",
                  variable_name))
}


