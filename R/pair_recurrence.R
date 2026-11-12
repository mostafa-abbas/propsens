#' Pair recurrence across models and across bootstrap resamples
#'
#' A descriptive index of how often the same untreated partner is chosen. Pair instability
#' (exploratory), reported by [propensity_sensitivity()], is 1 minus its mean across the kept
#' models over treated patients; it ranges from 0 to `1 - 1 / M` with `M` kept models. The
#' across-models and bootstrap versions are kept separate and are never mixed.
#'
#' **Across models.** For treated patient i, the modal partner is the untreated partner chosen
#' most often by the models in use; the recurrence is the share of those models choosing it (see
#' the denominator rule below). Floor `1 / M` (every model picks a different untreated partner),
#' ceiling 1.
#'
#' **Across bootstrap resamples** (when the match set has a bootstrap block: logistic regression
#' only, so it measures sampling instability). The modal partner is taken over the resamples.
#' The corrected recurrence divides the number of resamples in which i is matched to its modal
#' partner by the number of resamples in which **both** i and that partner are present; its
#' ceiling is 1. The uncorrected version divides by the number of resamples in which i is
#' present; its expected ceiling is about `1 - 1/e = 0.632`, because a given partner is missing
#' from about 37 % of bootstrap samples. Skipped resamples are ignored.
#'
#' Expect the across-models version near its floor (pair instability near its maximum) whenever
#' many untreated patients have similar scores, even when the models agree almost perfectly: the
#' index measures crowding more than model dependence. A different partner with the same
#' characteristics is harmless.
#'
#' @section Tie rule and denominator rule:
#' * **Ties.** When two or more untreated partners are chosen equally often, the modal partner is
#'   the one that comes first in the data (smallest row position). The same rule is used across
#'   models and across bootstrap resamples. With three models that all choose different
#'   untreated partners, the recurrence is 1/3 whichever partner is called modal.
#' * **Denominator across models.** The share is taken over the models in use that matched the
#'   patient. Without a caliper every model matches every treated patient, so the denominator
#'   is M, the number of models in use. With a caliper a model may leave a patient unmatched;
#'   that model is left out of the denominator, and a patient unmatched by every model gets
#'   `NA`. The reported floor is `1 / M`.
#' * **Denominator across bootstrap resamples.** Corrected: resamples in which both the patient
#'   and its modal partner are present. Uncorrected: resamples in which the patient is present.
#'   Skipped resamples are left out of both.
#'
#' @inheritSection rank_disagreement Quality gate
#' @param match_set A `match_set` from [match_each_model()].
#' @param use_gate `TRUE` (default) to use the models that passed the log-loss gate; `FALSE` to
#'   use every model. The bootstrap version does not depend on the gate.
#' @param models Optional character vector of model ids to use instead of a gate; it
#'   overrides `use_gate`. The wrapper [propensity_sensitivity()] passes the final
#'   quality-gate set (log-loss and balance) this way.
#' @return An object of class `"pair_recurrence"`: a list with `patients` (one row per patient:
#'   `row`, `id`, `treated`, `modal_control_models`, `pair_recurrence_models`,
#'   `modal_control_bootstrap`, `pair_recurrence_bootstrap`,
#'   `pair_recurrence_bootstrap_uncorrected`; `NA` for untreated patients), `summary` (means
#'   over treated,
#'   `floor_models`, `ceiling_models`, `ceiling_bootstrap`,
#'   `expected_ceiling_bootstrap_uncorrected`, `bootstrap_samples_used`), `models_used`,
#'   `quality_gate_left_fewer_than_two_models` and `settings`.
#' @examples
#' data("lalonde", package = "MatchIt")
#' model_set <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1, n_trees = 100,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
#' )
#' match_set <- match_each_model(model_set, lalonde, replace = TRUE, bootstrap_samples = 20)
#' recurrence <- pair_recurrence(match_set)
#' recurrence
#' @export
pair_recurrence <- function(match_set, use_gate = TRUE, models = NULL) {
  check_match_set(match_set)
  check_true_false(use_gate, "use_gate")
  n <- length(match_set$treatment)
  in_use <- models_in_use(match_set, use_gate, models)
  ids <- in_use$model_ids
  positions <- control_positions(match_set, ids)
  treated_rows <- which(match_set$treatment == 1)

  modal_models <- rep(NA_integer_, n)
  recurrence_models <- rep(NA_real_, n)
  for (i in treated_rows) {
    chosen <- positions[i, !is.na(positions[i, ])]
    if (length(chosen) == 0) next
    modal_models[i] <- modal_value(chosen)
    recurrence_models[i] <- mean(chosen == modal_models[i])
  }

  modal_bootstrap <- rep(NA_integer_, n)
  recurrence_bootstrap <- rep(NA_real_, n)
  recurrence_uncorrected <- rep(NA_real_, n)
  samples_used <- 0L
  if (!is.null(match_set$bootstrap)) {
    present <- match_set$bootstrap$present
    valid <- !is.na(present[1, ])
    samples_used <- sum(valid)
    present <- present[, valid, drop = FALSE]
    pairs <- matrix(match(match_set$bootstrap$pairs[, valid], match_set$row_ids), nrow = n)
    for (i in treated_rows) {
      chosen <- pairs[i, ]
      matched <- chosen[!is.na(chosen)]
      if (length(matched) == 0) next
      modal <- modal_value(matched)
      hits <- sum(chosen == modal, na.rm = TRUE)
      modal_bootstrap[i] <- modal
      recurrence_bootstrap[i] <- hits / sum(present[i, ] & present[modal, ])
      recurrence_uncorrected[i] <- hits / sum(present[i, ])
    }
  }

  treated <- match_set$treatment == 1
  mean_treated <- function(values) {
    if (all(is.na(values[treated]))) NA_real_ else mean(values[treated], na.rm = TRUE)
  }
  patients <- data.frame(
    row = seq_len(n), id = match_set$row_ids, treated = match_set$treatment,
    modal_control_models = match_set$row_ids[modal_models],
    pair_recurrence_models = recurrence_models,
    modal_control_bootstrap = match_set$row_ids[modal_bootstrap],
    pair_recurrence_bootstrap = recurrence_bootstrap,
    pair_recurrence_bootstrap_uncorrected = recurrence_uncorrected,
    stringsAsFactors = FALSE
  )
  structure(list(
    patients = patients,
    summary = list(mean_models_treated = mean_treated(recurrence_models),
                   floor_models = 1 / length(ids), ceiling_models = 1,
                   mean_bootstrap_treated = mean_treated(recurrence_bootstrap),
                   mean_bootstrap_uncorrected_treated = mean_treated(recurrence_uncorrected),
                   ceiling_bootstrap = 1,
                   expected_ceiling_bootstrap_uncorrected = 1 - exp(-1),
                   bootstrap_samples_used = samples_used),
    models_used = ids,
    quality_gate_left_fewer_than_two_models = in_use$quality_gate_left_fewer_than_two_models,
    settings = list(use_gate = use_gate, gate_applied = in_use$gate_applied,
                    models_supplied = in_use$models_supplied,
                    ties = "first control in the data")
  ), class = "pair_recurrence")
}

#' @rdname pair_recurrence
#' @param x A `pair_recurrence` object.
#' @param ... Not used.
#' @export
print.pair_recurrence <- function(x, ...) {
  print_diagnostic_header("Pair recurrence", x)
  s <- x$summary
  cat("Across models: mean over treated ", format_number(s$mean_models_treated),
      " (floor ", format_number(s$floor_models), ", ceiling 1)\n", sep = "")
  if (s$bootstrap_samples_used > 0) {
    cat("Across ", s$bootstrap_samples_used, " bootstrap resamples: corrected ",
        format_number(s$mean_bootstrap_treated), " (ceiling 1); uncorrected ",
        format_number(s$mean_bootstrap_uncorrected_treated), " (expected ceiling ",
        format_number(s$expected_ceiling_bootstrap_uncorrected), ")\n", sep = "")
  } else {
    cat("No bootstrap block in the match set.\n")
  }
  invisible(x)
}
