#' Repeat the diagnostic with several values of `seed`
#'
#' The fold split and every tree model of [propensity_sensitivity()] depend on `seed`, and on
#' small data the kept models and the summaries can change noticeably from one seed to the next.
#' `seed_sensitivity()` reruns [propensity_sensitivity()] once per value of `seeds` (here: values
#' of its `seed` argument) with every other argument unchanged, and returns one row per seed with
#' the models kept, whether the gate fell back, and the main summaries, plus their range. Report
#' the whole table.
#'
#' The R warnings of each run about the gate fallback and about a single kept algorithm are
#' not repeated; they are columns of the table (`fallback`, `algorithms_kept`).
#'
#' @param x A `matchit` object or a two-sided formula, as for [propensity_sensitivity()].
#' @param seeds Values of `seed` to run (default `1:3`). Each value is one full run.
#' @param ... Further arguments for [propensity_sensitivity()] (for example `data`,
#'   `algorithms`, `replace`, `caliper`, `cores`), except `seed`.
#' @param model_seeds The `seeds` argument of [propensity_sensitivity()] (the number of seeds
#'   per tree algorithm, as `1:k`); `NULL` (default) keeps its default `1:3`.
#' @param keep_results `TRUE` to keep every full result in `results` (large); default `FALSE`.
#' @return An object of class `"seed_sensitivity"`: a list with `table` (one row per seed:
#'   `seed`, `models_kept`, `models_fitted`, `algorithms_kept`, `kept_by_algorithm`,
#'   `fallback`, `rank_disagreement`, `rank_disagreement_equal_weight`, `pair_instability`,
#'   `bootstrap_rank_instability`), `range` (minimum and maximum of each summary over the seeds)
#'   and `results` (`NULL` unless `keep_results = TRUE`).
#' @examples
#' data("lalonde", package = "MatchIt")
#' by_seed <- seed_sensitivity(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75, seeds = 1:2,
#'   data = lalonde, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
#'   model_seeds = 1, n_trees = 50, bootstrap_samples = 2
#' )
#' by_seed
#' @export
seed_sensitivity <- function(x, seeds = 1:3, ..., model_seeds = NULL, keep_results = FALSE) {
  arguments <- list(...)
  # `seed = ` would be partially matched to `seeds`; it is refused, as is `seed` in `...`.
  if ("seed" %in% c(names(arguments), names(sys.call()))) {
    stop("Give the values of `seed` with `seeds`; do not pass `seed`.")
  }
  if (!is.numeric(seeds) || length(seeds) < 1 || anyNA(seeds) || any(seeds != round(seeds))) {
    stop("`seeds` must be a vector of integers (the values of `seed` to run).")
  }
  check_true_false(keep_results, "keep_results")
  if (!is.null(model_seeds)) arguments$seeds <- model_seeds
  # The two warnings below are columns of the table; any other warning reaches the user.
  table_warning <- function(w) {
    text <- conditionMessage(w)
    if (startsWith(text, "WARNING: fewer than two models passed") ||
          grepl("are random seeds of one algorithm", text, fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  }
  results <- lapply(seeds, function(value) {
    withCallingHandlers(do.call(propensity_sensitivity, c(list(x), arguments,
                                                          list(seed = value))),
                        warning = table_warning)
  })
  rows <- lapply(seq_along(seeds), function(i) {
    result <- results[[i]]
    table <- result$summary$table
    value <- function(name) table$with_gate[table$quantity == name]
    models <- result$models[!result$models$model_id %in% "reference", , drop = FALSE]
    kept <- models[models$model_id %in% result$summary$models_used, , drop = FALSE]
    counts <- table(factor(kept$algorithm, levels = unique(kept$algorithm)))
    data.frame(
      seed = seeds[i],
      models_kept = nrow(kept),
      models_fitted = nrow(models),
      algorithms_kept = length(counts),
      kept_by_algorithm = paste0(algorithm_display_names(names(counts)), " ",
                                 as.integer(counts), collapse = ", "),
      fallback = isTRUE(result$quality_gate_left_fewer_than_two_models),
      rank_disagreement = value("mean_rank_disagreement_treated"),
      rank_disagreement_equal_weight = value("mean_rank_disagreement_equal_weight_treated"),
      pair_instability = value("pair_instability_models_mean_treated"),
      bootstrap_rank_instability = value("bootstrap_rank_instability_treated"),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  measures <- c("rank_disagreement", "rank_disagreement_equal_weight", "pair_instability",
                "bootstrap_rank_instability")
  range_of <- function(values) {
    if (all(is.na(values))) c(NA_real_, NA_real_) else range(values, na.rm = TRUE)
  }
  ranges <- data.frame(summary = measures,
                       minimum = vapply(measures, function(m) range_of(table[[m]])[1],
                                        numeric(1)),
                       maximum = vapply(measures, function(m) range_of(table[[m]])[2],
                                        numeric(1)),
                       stringsAsFactors = FALSE)
  rownames(ranges) <- NULL
  structure(list(table = table, range = ranges,
                 results = if (keep_results) results else NULL),
            class = "seed_sensitivity")
}

#' @rdname seed_sensitivity
#' @param digits Significant digits for printed numbers.
#' @export
print.seed_sensitivity <- function(x, digits = 3, ...) {
  cat("Propensity model sensitivity over ", nrow(x$table), " values of `seed`\n", sep = "")
  table <- x$table
  for (column in c("rank_disagreement", "rank_disagreement_equal_weight", "pair_instability",
                   "bootstrap_rank_instability")) {
    table[[column]] <- signif(table[[column]], digits)
  }
  print(table, row.names = FALSE)
  cat("\nRange over the seeds:\n")
  ranges <- x$range
  ranges$minimum <- signif(ranges$minimum, digits)
  ranges$maximum <- signif(ranges$maximum, digits)
  print(ranges, row.names = FALSE)
  if (any(x$table$fallback) || any(x$table$algorithms_kept == 1)) {
    cat(wrap_text(paste("Some seeds fell back to every model or kept a single algorithm; their",
                        "values do not describe a set of similarly good models from several",
                        "algorithms.")), sep = "\n")
  }
  cat(wrap_text("Report every row; the values depend on the seed as much as shown here."),
      sep = "\n")
  invisible(x)
}
