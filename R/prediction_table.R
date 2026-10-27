##' Table of polygon predictions
##'
##' Formats the polygon summaries returned by [predict.RiskMap()] as a
##' presentation-ready table, with one row per polygon.
##'
##' @param object A `RiskMap_prediction` object from [predict.RiskMap()] with
##'   `boundaries`.
##' @param label Name of the column of the boundaries used to name the rows.
##'   By default, the first character or factor column, or the polygon number
##'   when there is none.
##' @param aggregate Aggregations to include, such as `"mean"` or `"sum"`. By
##'   default, all those calculated.
##' @param summaries Summaries to include, such as `"mean"`, `"exceedance"` or
##'   `"sd"`. By default, all those calculated.
##' @param quantity Predictive quantity to tabulate. Defaults to `"response"`.
##' @param sort_by Name of a summary, such as `"mean"`, used to sort the rows
##'   in decreasing order, using the first aggregation that has it. By default,
##'   rows keep the order of the boundaries.
##' @param digits Number of decimal places for values that are not counts.
##' @param col_names Optional column names, replacing the default headers.
##' @param ... Additional arguments passed to [knitr::kable()].
##'
##' @details
##' Each aggregation gives one or more columns. When the predictive mean and
##' both interval limits are available, they are combined in one column as
##' "mean (lower to upper)". Exceedance probabilities are omitted for sums of a
##' probability, because the thresholds refer to the probability itself.
##' Probabilities, including prevalence averaged over
##' a polygon and exceedance probabilities, are shown as percentages. Sums are
##' shown to three significant digits with thousands separators. The default
##' headers name the quantity as in the maps drawn by
##' [plot.RiskMap_prediction()].
##'
##' @return An object of class `knitr_kable`.
##' @seealso [predict.RiskMap()], [plot.RiskMap_prediction()]
##' @examples
##' data(italy_sim)
##' example_data <- italy_sim[1:50, ]
##' fit <- glgpm(y ~ gp(), data = example_data, family = "gaussian",
##'              messages = FALSE)
##' areas <- sf::st_sf(name = c("West", "East"),
##'                    geometry = sf::st_make_grid(example_data, n = c(2, 1)))
##' prediction <- predict(fit, newdata = example_data, boundaries = areas,
##'                       aggregate = "mean",
##'                       summaries = c("mean", "lower", "upper"),
##'                       messages = FALSE)
##' to_table(prediction)
##' @method to_table RiskMap_prediction
##' @export
to_table.RiskMap_prediction <- function(object, label = NULL,
                                        aggregate = NULL, summaries = NULL,
                                        quantity = "response", sort_by = NULL,
                                        digits = 1, col_names = NULL, ...) {
  x <- object
  if (is.null(x$boundaries) || is.null(x$boundary_summary)) {
    stop("to_table() needs polygon summaries; supply 'boundaries' to ",
         "predict()", call. = FALSE)
  }
  check_positive_integer(digits, "digits", allow_zero = TRUE)
  rows <- x$boundary_summary[x$boundary_summary$quantity == quantity, ,
                             drop = FALSE]
  if (!nrow(rows)) {
    stop("No polygon summaries for quantity '", quantity, "'", call. = FALSE)
  }
  aggregates <- unique(rows$aggregate)
  if (!is.null(aggregate)) {
    unknown <- setdiff(aggregate, aggregates)
    if (length(unknown)) {
      stop("Aggregation not calculated: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    aggregates <- aggregate
  }

  n <- nrow(x$boundaries)
  boundary_data <- st_drop_geometry(x$boundaries)
  if (is.null(label)) {
    text_columns <- names(boundary_data)[vapply(
      boundary_data, function(v) is.character(v) || is.factor(v), logical(1)
    )]
    row_labels <- if (length(text_columns)) {
      as.character(boundary_data[[text_columns[1]]])
    } else {
      as.character(seq_len(n))
    }
  } else {
    if (!is.character(label) || length(label) != 1L ||
        !label %in% names(boundary_data)) {
      stop("'label' must name a column of the boundaries", call. = FALSE)
    }
    row_labels <- as.character(boundary_data[[label]])
  }

  values_for <- function(aggregate_name, summary_name, threshold = NA) {
    keep <- rows$aggregate == aggregate_name & rows$summary == summary_name
    if (!is.na(threshold)) keep <- keep & rows$threshold %in% threshold
    out <- rep(NA_real_, n)
    out[rows$boundary_id[keep]] <- rows$value[keep]
    out
  }

  table <- data.frame(Area = row_labels, check.names = FALSE)
  sort_values <- NULL
  level <- paste0(format(100 * (x$level %||% 0.95)), "% interval")
  summary_order <- c("mean", "median", "sd", "cv", "interval_width",
                     "exceedance")
  for (aggregate_name in aggregates) {
    present <- unique(rows$summary[rows$aggregate == aggregate_name])
    if (!is.null(summaries)) present <- intersect(present, summaries)
    combined <- all(c("mean", "lower", "upper") %in% present)
    shown <- c(intersect(summary_order, present),
               if (!combined) intersect(c("lower", "upper"), present),
               setdiff(present, c(summary_order, "lower", "upper")))
    percent <- prediction_is_probability(x, quantity, "mean") &&
      aggregate_name %in% c("mean", "median", "min", "max")
    count <- identical(aggregate_name, "sum")
    format_value <- function(v, probability = percent) {
      out <- if (probability) {
        formatC(100 * v, format = "f", digits = as.integer(digits))
      } else if (count) {
        format(signif(v, 3), big.mark = ",", scientific = FALSE,
               trim = TRUE)
      } else {
        formatC(v, format = "f", digits = as.integer(digits))
      }
      out[is.na(v)] <- ""
      out
    }
    name <- prediction_table_quantity(x, quantity, aggregate_name)
    for (summary_name in shown) {
      if (summary_name == "exceedance") {
        # Thresholds are on the scale of the quantity, so they do not apply
        # to sums of a probability, such as numbers of people.
        if (count && prediction_is_probability(x, quantity, "mean")) next
        thresholds <- sort(unique(rows$threshold[
          rows$aggregate == aggregate_name & rows$summary == "exceedance"
        ]))
        for (threshold in thresholds) {
          shown_threshold <- if (percent) {
            paste0(format(100 * threshold, drop0trailing = TRUE), "%")
          } else {
            format(threshold)
          }
          header <- paste0("Probability ", tolower(name), " > ",
                           shown_threshold, " (%)")
          table[[header]] <- format_value(
            values_for(aggregate_name, "exceedance", threshold), TRUE
          )
        }
        next
      }
      v <- values_for(aggregate_name, summary_name)
      if (identical(summary_name, sort_by) && is.null(sort_values)) {
        sort_values <- v
      }
      unit <- if (percent && summary_name != "cv") "%" else NULL
      if (summary_name == "mean" && combined) {
        cells <- sprintf("%s (%s to %s)", format_value(v),
                         format_value(values_for(aggregate_name, "lower")),
                         format_value(values_for(aggregate_name, "upper")))
        header <- paste0(name, " (", paste(c(unit, level), collapse = ", "),
                         ")")
      } else {
        cells <- format_value(v, percent && summary_name != "cv")
        header <- paste0(prediction_table_summary(summary_name, name, level),
                         if (!is.null(unit)) " (%)")
      }
      table[[header]] <- cells
    }
  }

  if (!is.null(sort_by)) {
    if (is.null(sort_values)) {
      stop("'sort_by' must name a summary in the table, such as \"mean\"",
           call. = FALSE)
    }
    table <- table[order(-sort_values), , drop = FALSE]
  }
  dots <- list(...)
  if (!is.null(col_names)) {
    if (length(col_names) != ncol(table)) {
      stop("'col_names' must have ", ncol(table), " elements", call. = FALSE)
    }
    dots$col.names <- col_names
  }
  if (is.null(dots$row.names)) dots$row.names <- FALSE
  if (is.null(dots$align)) dots$align <- c("l", rep("r", ncol(table) - 1L))
  do.call(kable, c(list(x = table), dots))
}

#' Name of the areal quantity in table headers
#' @noRd
prediction_table_quantity <- function(x, quantity, aggregate) {
  label <- prediction_quantity_label(x, quantity)
  weighted <- !is.null(x$weights)
  name <- switch(
    aggregate,
    mean = paste(if (weighted) "Weighted mean" else "Mean", label),
    median = paste(if (weighted) "Weighted median" else "Median", label),
    sum = paste(if (weighted) "Weighted sum of" else "Sum of", label),
    min = paste("Minimum", label),
    max = paste("Maximum", label),
    paste(aggregate, label)
  )
  paste0(toupper(substring(name, 1L, 1L)), substring(name, 2L))
}

#' Header for one summary of an areal quantity
#' @noRd
prediction_table_summary <- function(summary, name, level) {
  lower_name <- tolower(substring(name, 1L, 1L))
  name_inside <- paste0(lower_name, substring(name, 2L))
  switch(
    summary,
    mean = name,
    median = paste("Median of", name_inside),
    lower = paste("Lower limit of", level, "for", name_inside),
    upper = paste("Upper limit of", level, "for", name_inside),
    sd = paste("SD of", name_inside),
    cv = paste("CV of", name_inside),
    interval_width = paste("Width of", level, "for", name_inside),
    paste(summary, "of", name_inside)
  )
}
