#' Summerise coverage
#'
#' This function creates a summary of list coverage for a given dataset.
#'
#' The dataset must include the list and cell indicator columns \code{L1},
#' \code{L2}, \code{Y11}, \code{Y10}, \code{Y01}, and \code{Y00}. In the
#' presence of overcoverage, it must also include the target inclusion indicator
#' column \code{in_target}.
#'
#' @param dataset A data.frame with cell and list inclusion indicators.
#' 
#' @return A named list with four elements: 
#'   \code{out_target_summary} (if dataset has list overcoverage);
#'   \code{in_target_coverage_summary};
#'   \code{in_target_cell_counts}; and
#'   \code{in_target_cell_proportions}.
#'
#' @examples
#' n_union  <- 1200
#' n_target <- 1000
#' 
#' L1 <- sample(c(TRUE, FALSE), n_target, prob = c(0.6, 0.4), replace = TRUE)
#' L2 <- sample(c(TRUE, FALSE), n_target, prob = c(0.8, 0.2), replace = TRUE)
#'
#' ## Dataset without overcoverage:
#' target_dataset <- data.frame(L1  = L1,
#'                              L2  = L2,
#'                              Y11 = L1 & L2,
#'                              Y10 = L1 & !L2,
#'                              Y01 = !L1 & L2,
#'                              Y00 = !L1 & !L2)
#' 
#' make_coverage_summary(target_dataset)
#' 
#' ## Dataset with overcoverage:
#' union_dataset <- rbind(target_dataset,
#'                        data.frame(L1  = rep(FALSE, n_union - n_target),
#'                                   L2  = TRUE,
#'                                   Y11 = NA,
#'                                   Y10 = NA,
#'                                   Y01 = NA,
#'                                   Y00 = NA))
#' 
#' union_dataset[["in_target"]] <- c(rep(TRUE, n_target), 
#'                                  rep(FALSE, n_union - n_target))
#'                                   
#' make_coverage_summary(union_dataset)
#'
#' @export
make_coverage_summary <- function(dataset) {
  ## Identify records in/outside target
  in_target  <- if ("in_target" %in% names(dataset)) 
                  dataset[["in_target"]] else 
                    rep(TRUE, nrow(dataset))
  
  indicators <- c("L1", "L2", "Y11", "Y10", "Y01", "Y00")
  dataset    <- dataset[indicators]
  n_target   <- sum(in_target)
  n_outside  <- sum(!in_target)

  ## List/cell counts in and outside target
  in_counts  <- c(colSums(dataset[in_target, ]), 
                  N_in_target = n_target)
  out_counts <- c(L1 = sum(dataset[!in_target, "L1"]),
                  L2 = sum(dataset[!in_target, "L2"]),
                  N_out_target = n_outside)
  
  ## Converting to factor to control sort order in xtabs()
  factors <- lapply(dataset[in_target, ],
                    function(x) 
                      factor(x, levels = c(TRUE, FALSE), labels = c("1", "0")))
  
  list(## outside target:
       out_target_summary = 
         data.frame(label      = names(out_counts),
                    count      = unname(out_counts),
                    proportion = format(if (n_outside > 0)
                                          round(out_counts / n_outside, 4) else 0,
                                        nsmall = 4),
                    row.names  = NULL),
       ## in target: summary
       in_target_coverage_summary = 
         data.frame(label      = names(in_counts),
                    count      = unname(in_counts),
                    proportion = format(round(in_counts / n_target, 4), 
                                        nsmall = 4),
                    row.names  = NULL),
       ## in target: cross-tabs
       in_target_cell_counts = 
         stats::xtabs(~ L1 + L2, data = data.frame(factors)),
       in_target_cell_proportions = 
         round(stats::xtabs(~ L1 + L2, data = data.frame(factors)) / n_target, 4)
  )
}

