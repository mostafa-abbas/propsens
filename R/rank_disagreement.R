#' Per-patient rank disagreement across propensity models
#'
#' For each patient, `d` is the standard deviation over
#' the models in use of the percentile rank `r_m(i) = rank(logit_m(i)) / n`. The decomposition
#' separates disagreement between algorithms from seed noise within an algorithm:
#' `d_between` is the standard deviation over algorithms of the per-algorithm mean rank, and
#' `d_within` is the square root of the mean over algorithms with at least two seeds of the
#' variance over seeds (`NA` when no algorithm has two seeds). A ratio `d_between / d_within`
#' near 1 means the algorithms disagree with themselves as much as with each other: noise, not
#' model dependence. `d_logit`, the standard deviation of the logits, is the secondary
#' logit-scale version; calibration differences inflate it without changing who is near whom.
#' Standard deviations use the n - 1 denominator.
#'
#' @section Quality gate:
#' With `use_gate = TRUE` only the models that passed the log-loss gate are used. If fewer than
#' two passed, all models are used and `quality_gate_left_fewer_than_two_models` is `TRUE`
#'. The flag reports the gate's outcome whatever `use_gate` is.
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
#'   `mean_d_all`, `median_between_within_ratio_treated`, `mean_d_logit_treated`), `models_used`,
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
    summary = list(mean_d_treated = mean(d[treated]), mean_d_all = mean(d),
                   median_between_within_ratio_treated = median_ratio,
                   mean_d_logit_treated = mean(d_logit[treated])),
    models_used = ids,
    quality_gate_left_fewer_than_two_models = in_use$quality_gate_left_fewer_than_two_models,
    settings = list(use_gate = use_gate, gate_applied = in_use$gate_applied,
                    models_supplied = in_use$models_supplied,
                    algorithms_with_seeds = with_seeds)
  ), class = "rank_disagreement")
}

#' @rdname rank_disagreement
#' @param x A `rank_disagreement` object.
#' @param ... Not used.
#' @export
print.rank_disagreement <- function(x, ...) {
  print_diagnostic_header("Rank disagreement", x)
  s <- x$summary
  cat("Mean d over treated: ", format_number(s$mean_d_treated),
      "; over all patients: ", format_number(s$mean_d_all), "\n", sep = "")
  cat("Median d_between / d_within over treated: ",
      format_number(s$median_between_within_ratio_treated),
      " (near 1 = seed noise, not model dependence)\n", sep = "")
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
