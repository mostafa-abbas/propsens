#' Per-patient rank disagreement across propensity models
#'
#' For each patient, `d` is the standard deviation over
#' the kept models of the percentile rank `r_m(i) = rank(logit_m(i)) / n` (0 to 1 scale).
#' Rank disagreement, the study mean of `d` over treated patients, is the primary summary; it is
#' a study-level summary, not a per-patient flag.
#'
#' Two exploratory companions, which are not a decomposition of `d`
#' (`d^2` is not `d_between^2 + d_within^2`):
#' `d_between` is the standard deviation over algorithms of each algorithm's mean rank, so that
#' every algorithm counts once however many seeds it has; its study mean over treated patients is
#' reported as "rank disagreement, equal weight per algorithm" (`NA` with one algorithm).
#' `d_within` is the square root of the mean, over algorithms with at least two seeds, of the
#' variance over seeds (`NA` when no algorithm has two seeds). The ratio `d_between / d_within` is
#' the disagreement between algorithms divided by the disagreement between seeds of the same
#' algorithm. With three seeds per algorithm, seed noise alone gives a ratio of about 0.4;
#' larger values point to disagreement between algorithms beyond seed noise.
#'
#' `d` is a standard deviation over models, so it depends on how many seeds each algorithm
#' contributes: with the same algorithms, more seeds of one algorithm pull `d` towards that
#' algorithm's seed noise. Values are comparable only between runs with the same model set
#' and the same kept algorithms. `d_logit`, the standard deviation of the logits, is a
#' logit-scale version; calibration differences inflate it without changing who is near whom.
#' Standard deviations use the n - 1 denominator.
#'
#' @section Quality gate:
#' "Kept models" in this help page are the models the function uses: the models given in
#' `models` (as [propensity_sensitivity()] passes the models kept by its quality gate), or with
#' `use_gate = TRUE` the models that passed the log-loss gate, or every model. If fewer than
#' two passed, all models are used and `quality_gate_left_fewer_than_two_models` is `TRUE`.
#' The flag reports the gate's outcome whatever `use_gate` is.
#'
#' @param model_set A `propensity_model_set` from [fit_propensity_models()].
#' @param use_gate `TRUE` (default) to use the models that passed the log-loss gate; `FALSE` to
#'   use every model.
#' @param models Optional character vector of model ids to use instead of a gate; it
#'   overrides `use_gate`. The wrapper [propensity_sensitivity()] passes the final
#'   quality-gate set (log-loss and balance) this way.
#' @return An object of class `"rank_disagreement"`: a list with `patients` (one row per patient:
#'   `row`, `treated`, `rank_disagreement`, `rank_disagreement_between`,
#'   `rank_disagreement_within`, `logit_disagreement`), `summary` (`mean_d_treated`,
#'   `mean_d_between_treated`, `mean_d_all`, `median_between_within_ratio_treated`,
#'   `mean_d_logit_treated`), `models_used`,
#'   `quality_gate_left_fewer_than_two_models` and `settings`.
#' @examples
#' data("lalonde", package = "MatchIt")
#' model_set <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1:2, n_trees = 100,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
#' )
#' disagreement <- rank_disagreement(model_set)
#' disagreement
#' head(disagreement$patients)
#' @export
rank_disagreement <- function(model_set, use_gate = TRUE, models = NULL) {
  if (!inherits(model_set, "propensity_model_set")) {
    stop("`model_set` must be a `propensity_model_set` from fit_propensity_models().")
  }
  check_true_false(use_gate, "use_gate")
  in_use <- models_in_use(model_set, use_gate, models)
  ids <- in_use$model_ids
  algorithms <- in_use$algorithms
  rank <- model_set$rank[, ids, drop = FALSE]
  logit <- model_set$logit[, ids, drop = FALSE]
  n <- nrow(rank)

  d <- row_sd(rank)
  d_logit <- row_sd(logit)
  algorithm_names <- unique(algorithms)
  algorithm_means <- matrix(vapply(algorithm_names, function(a) {
    rowMeans(rank[, algorithms == a, drop = FALSE])
  }, numeric(n)), nrow = n)
  d_between <- row_sd(algorithm_means)
  seed_counts <- table(algorithms)[algorithm_names]
  with_seeds <- algorithm_names[seed_counts >= 2]
  d_within <- if (length(with_seeds) == 0) {
    rep(NA_real_, n)
  } else {
    variances <- matrix(vapply(with_seeds, function(a) {
      row_sd(rank[, algorithms == a, drop = FALSE])^2
    }, numeric(n)), nrow = n)
    sqrt(rowMeans(variances))
  }

  treated <- model_set$treatment == 1
  ratio <- d_between / d_within
  median_ratio <- NA_real_
  if (!all(is.na(d_within))) median_ratio <- stats::median(ratio[treated], na.rm = TRUE)
  patients <- data.frame(
    row = seq_len(n), treated = as.integer(model_set$treatment),
    rank_disagreement = d, rank_disagreement_between = d_between,
    rank_disagreement_within = d_within, logit_disagreement = d_logit
  )
  structure(list(
    patients = patients,
    summary = list(mean_d_treated = mean(d[treated]),
                   mean_d_between_treated = mean(d_between[treated]),
                   mean_d_all = mean(d),
                   median_between_within_ratio_treated = median_ratio,
                   mean_d_logit_treated = mean(d_logit[treated])),
    models_used = ids,
    quality_gate_left_fewer_than_two_models = in_use$quality_gate_left_fewer_than_two_models,
    settings = list(use_gate = use_gate, gate_applied = in_use$gate_applied,
                    models_supplied = in_use$models_supplied,
                    algorithms_with_seeds = with_seeds,
                    seeds_per_algorithm = c(seed_counts))
  ), class = "rank_disagreement")
}

#' @rdname rank_disagreement
#' @param x A `rank_disagreement` object.
#' @param ... Not used.
#' @export
print.rank_disagreement <- function(x, ...) {
  print_diagnostic_header("Rank disagreement", x)
  s <- x$summary
  cat("Mean d over treated (primary): ", format_number(s$mean_d_treated),
      "; over all patients: ", format_number(s$mean_d_all), "\n", sep = "")
  cat("Equal weight per algorithm, mean d_between over treated (exploratory): ",
      format_number(s$mean_d_between_treated), "\n", sep = "")
  counts <- x$settings$seeds_per_algorithm
  tree_counts <- counts[!names(counts) %in% single_fit_algorithms]
  reference_applies <- !is.na(s$median_between_within_ratio_treated) &&
    length(tree_counts) > 0 && all(tree_counts == 3)
  cat(wrap_text(paste0("Median d_between / d_within over treated (exploratory): ",
                       format_number(s$median_between_within_ratio_treated),
                       if (reference_applies) {
                         paste0(". With three seeds per algorithm, seed noise alone gives about ",
                                "0.4; larger values point to disagreement between algorithms.")
                       })), sep = "\n")
  cat("Mean d_logit over treated: ", format_number(s$mean_d_logit_treated), "\n", sep = "")
  invisible(x)
}

print_diagnostic_header <- function(title, x) {
  cat(title, ": ", length(x$models_used), " model(s) in use", sep = "")
  if (isTRUE(x$settings$models_supplied)) {
    cat(" (models chosen by the caller)")
  } else if (isTRUE(x$quality_gate_left_fewer_than_two_models)) {
    cat(" (quality gate left fewer than two models: all models used)")
  } else if (isTRUE(x$settings$gate_applied)) {
    cat(" (after the log-loss gate)")
  } else {
    cat(" (gate not applied)")
  }
  cat("\n")
}

format_number <- function(value, digits = 4) {
  if (is.null(value) || length(value) == 0 || is.na(value)) return("NA")
  format(signif(value, digits))
}

# A tuning constant such as 0.1 printed with two decimals ("0.10"), as in the documentation.
format_setting <- function(value) {
  formatC(value, format = "f", digits = max(2, nchar(sub("^[^.]*\\.?", "", format(value)))))
}

# Lines wrapped to the console width; continuation lines are indented by `exdent` spaces.
wrap_text <- function(text, indent = 0, exdent = indent + 3, width = getOption("width", 80)) {
  strwrap(text, width = max(40, width - 1), indent = indent, exdent = exdent)
}
