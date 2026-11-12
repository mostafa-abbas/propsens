#' Methods and Results sentences for a propensity model sensitivity result
#'
#' Builds the sentences to report a [propensity_sensitivity()] result, from the result itself:
#' the model set, the quality gate and the models it kept (by algorithm), the matching settings,
#' the share of treated patients kept, the seed, and the main summaries. After a fallback (fewer
#' than two models passed the gate) the text says that all models were used and describes them
#' as all models. Edit the sentences to fit your paper; call the gate criteria pre-specified only
#' if you fixed them before running the check.
#'
#' @param x A `propensity_sensitivity` object.
#' @param digits Decimal places for the rank summaries.
#' @return A named character vector with `methods` and `results`, one paragraph each.
#' @examples
#' data("lalonde", package = "MatchIt")
#' sensitivity <- propensity_sensitivity(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
#'   seeds = 1, n_trees = 50, bootstrap_samples = 5, replace = TRUE
#' )
#' cat(methods_text(sensitivity), sep = "\n\n")
#' @export
methods_text <- function(x, digits = 3) {
  if (!inherits(x, "propensity_sensitivity")) {
    stop("`x` must be a propensity_sensitivity object.")
  }
  check_single_integer(digits, "digits", minimum = 1)
  settings <- x$settings
  models <- x$models[!x$models$is_reference, , drop = FALSE]
  kept <- models[models$in_use, , drop = FALSE]
  fallback <- isTRUE(x$quality_gate_left_fewer_than_two_models)
  table <- x$summary$table
  value <- function(name) table$with_gate[table$quantity == name]

  fitted_algorithms <- unique(models$algorithm)
  seed_text <- if (any(!fitted_algorithms %in% single_fit_algorithms)) {
    paste0(length(settings$seeds), if (length(settings$seeds) == 1) " seed" else " seeds",
           " for each tree method; ")
  } else {
    ""
  }
  balance_rule <- switch(
    settings$balance_rule,
    relative = sprintf("within %s of the best model's and at most %s",
                       format_setting(settings$balance_margin),
                       format_setting(settings$balance_cap)),
    relative_only = sprintf("within %s of the best model's",
                            format_setting(settings$balance_margin)),
    absolute = sprintf("at most %s", format_setting(settings$balance_threshold))
  )
  caliper <- if (is.null(settings$caliper_sd)) "" else
    sprintf(" with a caliper of %s standard deviations of the logit",
            format(signif(settings$caliper_sd, 2)))
  share <- range(kept$treated_kept)
  share_text <- if (all(share >= 1)) "all treated patients were matched" else
    sprintf("%.0f to %.0f %% of treated patients were matched", 100 * share[1], 100 * share[2])
  counts <- table(factor(kept$algorithm, levels = unique(kept$algorithm)))
  gate_result <- if (fallback) {
    sprintf(paste("fewer than two of the %d models met these criteria, so all models were",
                  "used"), nrow(models))
  } else {
    sprintf("%d of %d models were kept (%s)", nrow(kept), nrow(models),
            paste0(algorithm_display_names(names(counts)), ": ", as.integer(counts),
                   collapse = "; "))
  }
  which_models <- if (fallback) "all models" else "the kept models"

  methods <- paste(
    sprintf(paste("We described the sensitivity of the matched design to the propensity model",
                  "with the R package propsens (version %s)."),
            as.character(utils::packageVersion("propsens"))),
    sprintf(paste("Propensity models from %s (%s%d models; %d-fold cross-fitted scores; seed",
                  "%d) were each used for 1:1 nearest-neighbour matching %s replacement%s."),
            english_list(algorithm_display_names(fitted_algorithms)), seed_text, nrow(models),
            settings$folds, settings$seed, if (settings$replace) "with" else "without",
            caliper),
    sprintf(paste("Models were kept if their cross-fitted log-loss was within %s %% of the best",
                  "model's and their largest post-match standardized mean difference on main",
                  "terms was %s; %s, and %s."),
            format(100 * settings$quality_tolerance), balance_rule, gate_result, share_text),
    sprintf(paste("For %s we report the mean rank disagreement over treated patients, the",
                  "bootstrap rank instability of the logistic model, and the range of the",
                  "largest standardized mean difference on main terms."), which_models),
    paste("No outcome data were used. These summaries are descriptive; they have not been",
          "validated as predictors of bias.")
  )

  n_used <- length(x$summary$models_used)
  ratio <- value("median_between_within_ratio_treated")
  in_use <- models[models$model_id %in% x$summary$models_used, , drop = FALSE]
  number <- function(v) {
    if (is.na(v)) "not available" else formatC(v, format = "f", digits = digits)
  }
  sentences <- c(
    sprintf(paste("Mean rank disagreement over treated patients was %s across %s; the",
                  "bootstrap rank instability of the logistic model was %s."),
            number(value("mean_rank_disagreement_treated")), which_models,
            number(value("bootstrap_rank_instability_treated"))),
    if (n_used >= 2) {
      sprintf(paste("As exploratory summaries, pair instability was %.2f (of a possible %.2f",
                    "with %d models) and the median between/within-algorithm ratio was %s."),
              value("pair_instability_models_mean_treated"), 1 - 1 / n_used, n_used,
              if (is.na(ratio)) "not available" else sprintf("%.1f", ratio))
    },
    sprintf(paste("The largest standardized mean difference on main terms ranged from %.2f to",
                  "%.2f across %s."),
            min(in_use$smd_max_main), max(in_use$smd_max_main), which_models),
    if (isTRUE(x$summary$kept_models_one_algorithm)) {
      paste("The models used were seeds of a single algorithm, so these summaries describe",
            "seed noise within that algorithm, not disagreement between algorithms.")
    }
  )
  c(methods = methods, results = paste(sentences, collapse = " "))
}

# "a, b and c".
english_list <- function(words) {
  if (length(words) <= 1) return(paste(words))
  paste(paste(words[-length(words)], collapse = ", "), "and", words[length(words)])
}
