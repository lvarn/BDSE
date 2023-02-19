#' Make distribution over unique values
#'
#' Aggregates count columns over unique values of the specified grouping
#' columns.
#' 
#' When the count column is 1 for all rows, this is the frequency distribution
#' over the unique values of the grouping variables.
#' 
#' @param df A data.frame or data.table.
#' @param contract_cols A character vector giving columns to aggregate.
#' @param complete_cols A character vector giving the columns in \code{df} to
#'   keep which must include columns specified in \code{contract_cols}; defaults
#'   to \code{names(df)}.
#' @param as_dataframe A logical scalar, specifying whether the output should be
#'   coerced to a data.frame (defaults to \code{TRUE}). If \code{FALSE}, the
#'   function returns a data.table.
#'
#' @return A data.frame or data.table of aggragted counts.
#'
#' @import data.table
#'
#' @examples
#' n   <- 20
#' counts <- sample(0:1, n, replace = TRUE)
#' df  <- data.frame(x1 = sample(1:2, n, replace = TRUE),
#'                   x2 = letters[1:2],
#'                   y1 = counts,
#'                   y2 = 1 - counts)
#'
#' distribute(df, c("y1", "y2"))
#' distribute(df, c("y1"), c("x2", "y1"))
#'
#' @export
distribute <- function(df,
                       contract_cols,
                       complete_cols = names(df),
                       as_dataframe  = TRUE) {
  df <- as.data.table(df)[, complete_cols, with = FALSE]
  
  if (!all(stats::complete.cases(df)))
    message("Some columns include NAs")

  df <- df[
            stats::complete.cases(df), 
            lapply(.SD, sum, na.rm = TRUE),
            by = setdiff(names(df), contract_cols),
            .SDcols = contract_cols
          ]
  
  if (as_dataframe) df <- as.data.frame(df)
  
  df
}


