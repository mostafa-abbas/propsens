#' Per-patient control spread across propensity models
#'
#' For each treated patient, collect the distinct
#' untreated partners matched to it by the models in use and compute `s`, the mean pairwise
#' Mahalanobis distance between their covariate vectors, with the covariance matrix of all n
#' patients. `s = 0` when every model picks the same untreated partner, and `NA` when the
#' patient is unmatched in
#' every model. `s_scaled = s / D`, where `D` is the mean over matched treated patients of the
#' Mahalanobis distance between a treated patient and its logistic-regression match (the first
#' model in use when the match set has no logistic regression model): 1 means the candidate
#' partners differ from each other as much as a typical match differs from its treated patient.
#' The oracle entry and the bootstrap block are never used.
#'
#' @section Covariates:
#' Factors and logical columns are expanded to dummies (treatment contrasts). Character columns
#' are an error: convert them to factors. Constant columns and linearly dependent columns (found
#' with a pivoted QR decomposition of the standardized design) are dropped with a message. The
#' Mahalanobis distance is `sqrt((a - b)' S^-1 (a - b))`, the square root of what
#' [stats::mahalanobis()] returns.
#'
#' @inheritSection rank_disagreement Quality gate
#' @param match_set A `match_set` from [match_each_model()].
#' @param covariates A data frame with one row per patient, in the data's row order, holding the
#'   covariates that define the distance (usually those of the propensity formula).
#' @param use_gate `TRUE` (default) to use the models that passed the log-loss gate; `FALSE` to
#'   use every model.
#' @param models Optional character vector of model ids to use instead of a gate; it
#'   overrides `use_gate`. The wrapper [propensity_sensitivity()] passes the final
#'   quality-gate set (log-loss and balance) this way.
#' @return An object of class `"control_spread"`: a list with `patients` (one row per patient:
#'   `row`, `id`, `treated`, `distinct_controls`, `control_spread`, `control_spread_scaled`;
#'   `NA` for untreated patients), `summary` (`mean_scaled_treated`, `p90_scaled_treated`,
#'   `scale_distance`, `scale_model`), `models_used`,
#'   `quality_gate_left_fewer_than_two_models` and `settings` (with the dropped columns).
#' @examples
#' data("lalonde", package = "MatchIt")
#' model_set <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1, n_trees = 100,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
#' )
#' match_set <- match_each_model(model_set, lalonde, replace = TRUE)
#' covariates <- lalonde[, c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
#' spread <- control_spread(match_set, covariates)
#' spread
#' head(spread$patients[spread$patients$treated == 1, ])
#' @export
control_spread <- function(match_set, covariates, use_gate = TRUE, models = NULL) {
  check_match_set(match_set)
  check_true_false(use_gate, "use_gate")
  n <- length(match_set$treatment)
  prepared <- prepare_covariate_matrix(covariates, n, "control_spread")
  whitened <- whiten_rows(prepared$design)
  in_use <- models_in_use(match_set, use_gate, models)
  ids <- in_use$model_ids
  positions <- control_positions(match_set, ids)
  treated_rows <- which(match_set$treatment == 1)

  distinct <- rep(NA_integer_, n)
  spread <- rep(NA_real_, n)
  for (i in treated_rows) {
    controls <- unique(positions[i, !is.na(positions[i, ])])
    distinct[i] <- length(controls)
    if (length(controls) == 1) spread[i] <- 0
    if (length(controls) > 1) spread[i] <- mean(stats::dist(whitened[controls, , drop = FALSE]))
  }

  logistic_ids <- match_set$models$model_id[match_set$models$algorithm == "logistic_regression"]
  scale_model <- if (length(logistic_ids) > 0) logistic_ids[1] else ids[1]
  scale_pairs <- match_set$matches[[scale_model]]
  treated_position <- match(scale_pairs$treated_id, match_set$row_ids)
  control_position <- match(scale_pairs$control_id, match_set$row_ids)
  pair_distances <- sqrt(rowSums((whitened[treated_position, , drop = FALSE] -
                                    whitened[control_position, , drop = FALSE])^2))
  scale_distance <- mean(pair_distances)
  scaled <- spread / scale_distance

  treated <- match_set$treatment == 1
  patients <- data.frame(
    row = seq_len(n), id = match_set$row_ids, treated = match_set$treatment,
    distinct_controls = distinct, control_spread = spread, control_spread_scaled = scaled,
    stringsAsFactors = FALSE
  )
  structure(list(
    patients = patients,
    summary = list(mean_scaled_treated = mean(scaled[treated], na.rm = TRUE),
                   p90_scaled_treated = unname(stats::quantile(scaled[treated], 0.9,
                                                               na.rm = TRUE)),
                   scale_distance = scale_distance, scale_model = scale_model,
                   treated_unmatched_in_every_model = sum(treated & is.na(spread))),
    models_used = ids,
    quality_gate_left_fewer_than_two_models = in_use$quality_gate_left_fewer_than_two_models,
    settings = list(use_gate = use_gate, gate_applied = in_use$gate_applied,
                    models_supplied = in_use$models_supplied,
                    covariate_columns = colnames(prepared$design),
                    dropped_constant = prepared$dropped_constant,
                    dropped_dependent = prepared$dropped_dependent)
  ), class = "control_spread")
}

#' @rdname control_spread
#' @param x A `control_spread` object.
#' @param ... Not used.
#' @export
print.control_spread <- function(x, ...) {
  print_diagnostic_header("Control spread", x)
  s <- x$summary
  cat("Scaled control spread over treated: mean ", format_number(s$mean_scaled_treated),
      ", 90th percentile ", format_number(s$p90_scaled_treated), "\n", sep = "")
  cat("Scale: mean Mahalanobis distance treated to match under ", s$scale_model, " = ",
      format_number(s$scale_distance), "\n", sep = "")
  invisible(x)
}
