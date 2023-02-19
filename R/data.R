#' Synthetic base dataset
#'
#' This is an extension of \code{base_dataset_raw}. It includes individuals with
#' ages below 15, single-year ages for al records and ethnicity indicators.
#' 
#' @format A data.frame with 37108 rows and 20 variables:
#' \describe{
#'   \item{ethnicity}{...}
#'   \item{region}{...}
#'   \item{sex}{...}
#'   \item{age_grp}{...}
#'   \item{qualification}{...}
#'   \item{occupation}{...}
#'   \item{hours}{...}
#'   \item{income}{...}
#'   \item{hours_cen_mean}{...}
#'   \item{hours_cen_40hrs}{...}
#'   \item{hours_grp}{...}
#'   \item{age}{...}
#'   \item{region_fact}{...}
#'   \item{ethnicity1}{...}
#'   \item{ethnicity2}{...}
#'   \item{ethnicity3}{...}
#'   \item{ethnicity4}{...}
#'   \item{ethnicity5}{...}
#'   \item{ethnicity6}{...}
#'   \item{ethnicity9}{...}
#' }
"base_dataset"


#' Base dataset (From StatsNZ Website)
#'
#' A dataset containing covariate values, which can be used to generate a 
#' realistic synthetic target population dataset.
#'
#' @format A data.frame with 29471 rows and 8 variables:
#' \describe{
#'   \item{ethnicity}{...}
#'   \item{region}{...}
#'   \item{sex}{...} 
#'   \item{age_grp}{...}
#'   \item{qualification}{...}
#'   \item{occupation}{...}
#'   \item{hours}{...}
#'   \item{income}{...}
#' }
"base_dataset_raw"


#' Fertility data
#'
#' A dataset containing the fertility rates for women of childbearing age by age
#' group.
#'
#' @format A data.frame with 11 rows and 3 variables:
#' \describe{
#'   \item{mothers_age_grp}{...}
#'   \item{average_births}{...}
#'   \item{fert_rate}{...}
#' }
"fertility"

