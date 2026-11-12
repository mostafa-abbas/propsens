#' Propensity model sensitivity of a matched design
#'
#' The one-call diagnostic. It refits the propensity model with
#' every algorithm of the model set, re-matches with the analyst's matching options, applies the
#' quality gate, and reports rank disagreement (the primary summary), pair instability (the second
#' summary of the same gated set: 1 minus the pair recurrence across models; exploratory), the
#' bootstrap rank instability of the logistic model, control spread and pair recurrence, and
#' study-level balance spread. No outcome is used.
#' Rank disagreement is a study-level diagnostic; read it for the study as a whole, not per
#' patient. The diagnostics are evaluated in a manuscript in preparation.
#' When fewer than two models pass the quality gate, `print()` and `summary()` say so first: the
#' diagnostic is then outside the setting where it was evaluated.
#'
#' @section Steps:
#' 1. [fit_propensity_models()] with the formula (cross-fitted scores by default).
#' 2. [match_each_model()] with the matching options and the logistic bootstrap block.
#' 3. [balance_spread()] for every model.
#' 4. Final quality gate: `passed_quality_gate = passed_logloss_gate & smd_max_main <=
#'    balance_limit`. With the default `balance_rule = "relative"` the limit is
#'    the smallest max |SMD| on main terms among the refitted models plus `balance_margin`
#'    (0.10), and never above `balance_cap` (0.25): a model is excluded only when its balance is
#'    clearly worse than the best model's. With `balance_rule = "relative_only"` the limit
#'    is the best model's value plus `balance_margin`, without the cap; use it as a sensitivity
#'    analysis when the cap excludes every model (many covariates). With `balance_rule =
#'    "absolute"` the limit is `balance_threshold` (0.10). If fewer than two models pass,
#'    every model is used and
#'    `quality_gate_left_fewer_than_two_models` is `TRUE` (applied to the final set).
#' 5. [rank_disagreement()], [control_spread()], [pair_recurrence()] and the balance spread on
#'    the final set (`summary` column `with_gate`) and on every model (`without_gate`).
#'
#' @section Reading a matchit object:
#' From a `matchit` object the function reads the formula, the method, `ratio`, `replace`,
#' `m.order` (from the call; MatchIt's default for the ATT with a propensity score is
#' `"largest"`) and the caliper. Only 1:1 nearest-neighbour matching on a propensity score for
#' the ATT is supported, with or without replacement and with or without a caliper on the
#' score; anything else stops with a clear error. MatchIt applies a caliper in standard
#' deviations of the analyst's distance (over all patients); the same number of standard
#' deviations is applied to the logit of every model here. When `data` is `NULL` the data is
#' found by evaluating the object's `data` argument in the calling frame.
#'
#' @section Reference model:
#' The analyst's own score (`x$distance`) is kept as the reference model, an extra row of
#' `models` with `model_id = "reference"`. Its ranks come from `x$distance` (ranks do not depend
#' on the scale) and its balance statistics from `x$weights`. Its log-loss, Brier score and
#' c-statistic are `NA`: it is an in-sample score from the user's object and is not comparable
#' with the cross-fitted scores. Its quality gate is judged on balance only; when it fails, the
#' row is flagged "reference model outside the equally-reasonable set". It is never one of the
#' models in use. `summary` reports its mean absolute rank difference from the models in use,
#' over treated patients: the mean over treated `i` and models `m` of `|r_ref(i) - r_m(i)|`.
#'
#' @param x A `matchit` object (from `MatchIt::matchit()`) or a two-sided formula such as
#'   `treat ~ x1 + x2`. Only baseline covariates belong in the formula; exclude outcomes,
#'   identifiers and post-treatment variables (for example, `treat ~ . - death - id`).
#' @param data The data frame. Required with a formula; with a `matchit` object it is looked up
#'   from the object's call when `NULL`.
#' @param algorithms The propensity model set. The default is the set evaluated in a
#'   manuscript in preparation:
#'   logistic regression, BART and the GRF probability forest (with the default `seeds`, three
#'   seeds for each tree algorithm: seven models; packages 'dbarts' and 'grf'). A faster set
#'   is available with `algorithms = c("logistic_regression", "random_forest",
#'   "gradient_boosting")`; it is not part of that evaluation. See
#'   [fit_propensity_models()] for every option and the fixed settings.
#' @param seeds,folds,cross_fit,n_trees,quality_tolerance Passed to
#'   [fit_propensity_models()]. `n_trees` applies to gradient boosting only.
#' @param balance_rule The balance part of the quality gate: `"relative"` (default),
#'   `"relative_only"` (the relative rule without the cap) or `"absolute"`. See
#'   Steps.
#' @param balance_margin,balance_cap Relative rule: a model fails if its post-match maximum
#'   absolute standardized mean difference on main terms is more than `balance_margin` above
#'   the best model's (default 0.10) or above `balance_cap` (default 0.25; not used by
#'   `"relative_only"`).
#' @param balance_threshold Absolute rule: largest allowed post-match maximum absolute
#'   standardized mean difference on main terms (default 0.10). Used only with
#'   `balance_rule = "absolute"`.
#' @param bootstrap_samples Number of logistic bootstrap refits, for the bootstrap rank
#'   instability and the bootstrap pair recurrence (default 20; 0 for none).
#' @param covariates Data frame of covariates for the control spread and the balance
#'   statistics. Its rows must be the same patients in the same order as `data`; when both
#'   have row names other than the default and they differ, the function stops. `NULL`
#'   (default) uses the source columns of the formula's right-hand side, from `data`, resolved
#'   with R's formula rules: `.` means every other column and exclusions such as `- death`
#'   are respected; `I(age^2)` gives `age`, `age:diabetes` gives `age` and `diabetes`.
#' @param replace,m_order,caliper Matching options for the formula interface (without replacement,
#'   order `"largest"`, no caliper by default; the caliper is in standard deviations of each
#'   model's logit). With a `matchit` object they are read from the object and must not be
#'   given.
#' @param cores Number of cores (default 1).
#' @param seed Integer seed from which every other seed is derived.
#' @return An object of class `"propensity_sensitivity"`: a list with `patients` (one row per
#'   patient: `id`, `treated`, `rank_disagreement`, `rank_disagreement_between`,
#'   `rank_disagreement_within`, `logit_disagreement`, `control_spread`,
#'   `control_spread_scaled`, `pair_recurrence_models`, `pair_recurrence_bootstrap`,
#'   `outside_support`, `logistic_score`; the treated-only measures are `NA` for untreated
#'   patients),
#'   `models` (one row per model plus the reference row), `summary` (a table of the study-level
#'   numbers with and without the gate, plus flags), `quality_gate_left_fewer_than_two_models`,
#'   `balance_spread` (spread tables with and without the gate), `bootstrap_rank_instability`
#'   (see [bootstrap_rank_instability()]; `NULL` without bootstrap refits; its study-level value
#'   is also the row `bootstrap_rank_instability_treated` of the summary table; the summary
#'   table also has `pair_instability_models_mean_treated`, 1 minus
#'   `pair_recurrence_models_mean_treated`), `settings` and
#'   `components`
#'   (the model set and the match set). `components$model_set` is the fitted
#'   [fit_propensity_models()] object, with every model's cross-fitted `probability`, `logit`
#'   and `rank`, so that its scores can be reused (for example in another estimator) without
#'   refitting; `summary$models_used` names the models in use.
#' @examples
#' data("lalonde", package = "MatchIt")
#' matched <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
#'                             data = lalonde, distance = "glm", link = "linear.logit",
#'                             replace = TRUE)
#' # A fast model set for this example; the default set is slower (see below).
#' sensitivity <- propensity_sensitivity(
#'   matched, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
#'   seeds = 1, n_trees = 100, bootstrap_samples = 5
#' )
#' sensitivity
#' summary(sensitivity)
#' head(as.data.frame(sensitivity))
#' \donttest{
#' # The defaults: logistic regression, BART x 3, GRF probability forest x 3 (about a minute).
#' propensity_sensitivity(matched)
#' }
#' @export
propensity_sensitivity <- function(x, data = NULL,
                                   algorithms = c("logistic_regression", "bart",
                                                  "grf_probability"),
                                   seeds = 1:3, folds = 5, cross_fit = TRUE, n_trees = NULL,
                                   quality_tolerance = 0.05,
                                   balance_rule = c("relative", "absolute", "relative_only"),
                                   balance_margin = 0.10, balance_cap = 0.25,
                                   balance_threshold = 0.10, bootstrap_samples = 20,
                                   covariates = NULL, replace = FALSE,
                                   m_order = "largest", caliper = NULL, cores = 1, seed = 1) {
  algorithms <- match.arg(algorithms, propensity_algorithms, several.ok = TRUE)
  balance_rule <- match.arg(balance_rule)
  for (name in c("balance_margin", "balance_cap", "balance_threshold")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1 || is.na(value) || value < 0) {
      stop("`", name, "` must be a single number of at least 0.")
    }
  }
  rule <- list(rule = balance_rule, margin = balance_margin, cap = balance_cap,
               threshold = balance_threshold)
  caller <- parent.frame()
  if (inherits(x, "matchit")) {
    given <- c("replace", "m_order", "caliper")[c(!missing(replace), !missing(m_order),
                                                  !missing(caliper))]
    if (length(given) > 0) {
      stop("The matching options are read from the matchit object; do not pass ",
           paste0("`", given, "`", collapse = ", "), ".")
    }
    if (is.null(data)) {
      data <- tryCatch(eval(x$call$data, caller), error = function(e) NULL)
      if (!is.data.frame(data)) {
        stop("Could not find the data of the matchit object; pass the data with `data =`.")
      }
    }
    design <- read_matchit_design(x, data, caller)
  } else if (inherits(x, "formula")) {
    if (is.null(data)) stop("With a formula, pass the data with `data =`.")
    design <- list(formula = x, replace = replace, m_order = m_order, caliper = caliper,
                   reference = NULL, source = "formula", matchit_call = NULL)
  } else {
    stop("`x` must be a matchit object (from MatchIt::matchit()) or a formula.")
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  formula <- design$formula
  covariates <- select_covariates(covariates, formula, data)

  model_set <- fit_propensity_models(formula, data, algorithms = algorithms, seeds = seeds,
                                     folds = folds, cross_fit = cross_fit, n_trees = n_trees,
                                     quality_tolerance = quality_tolerance, cores = cores,
                                     seed = seed)
  match_set <- match_each_model(model_set, data, replace = design$replace, ratio = 1,
                                m_order = design$m_order, caliper = design$caliper,
                                bootstrap_samples = bootstrap_samples, cores = cores,
                                seed = seed)
  treatment <- model_set$treatment
  instability <- if (bootstrap_samples >= 2) bootstrap_rank_instability(match_set)
  instability_value <- if (is.null(instability)) NA_real_ else instability$summary$mean_treated
  balance_all <- balance_spread(match_set, treatment, covariates, use_gate = FALSE)

  gate <- apply_final_quality_gate(model_set$models, balance_all$models, rule)
  models <- gate$models
  models$treated_kept <- unname(match_set$treated_kept[models$model_id])
  all_ids <- models$model_id
  final_ids <- gate$final_ids

  with_gate <- run_set_diagnostics(model_set, match_set, covariates, final_ids)
  without_gate <- suppressMessages(run_set_diagnostics(model_set, match_set, covariates,
                                                       all_ids))
  spread_with_gate <- balance_spread_table(balance_all$models, final_ids)
  spread_without_gate <- balance_spread_table(balance_all$models, all_ids)

  reference_difference <- c(NA_real_, NA_real_)
  if (!is.null(design$reference)) {
    reference_row <- reference_model_row(design$reference, covariates, treatment,
                                         gate$balance_limit, names(models))
    models <- rbind(models, reference_row)
    treated <- treatment == 1
    reference_difference <- vapply(list(final_ids, all_ids), function(ids) {
      mean(abs(model_set$rank[treated, ids, drop = FALSE] - design$reference$rank[treated]))
    }, numeric(1))
  }
  rownames(models) <- NULL

  summary_with_gate <- set_summary(with_gate, spread_with_gate, models, final_ids)
  summary_without_gate <- set_summary(without_gate, spread_without_gate, models, all_ids)
  summary_table <- data.frame(
    quantity = c("number_of_models_used", names(summary_with_gate),
                 "bootstrap_rank_instability_treated",
                 "reference_mean_absolute_rank_difference_treated"),
    with_gate = c(length(final_ids), unlist(summary_with_gate), instability_value,
                  reference_difference[1]),
    without_gate = c(length(all_ids), unlist(summary_without_gate), instability_value,
                     reference_difference[2]),
    stringsAsFactors = FALSE
  )
  rownames(summary_table) <- NULL

  reference_outside <- NA
  if (!is.null(design$reference)) {
    reference_outside <- !isTRUE(models$passed_quality_gate[models$model_id == "reference"])
  }
  structure(list(
    patients = patient_table(with_gate, model_set, match_set),
    models = models,
    summary = list(
      table = summary_table,
      quality_gate_left_fewer_than_two_models = gate$fewer_than_two,
      reference_model_outside_equally_reasonable_set = reference_outside,
      models_used = final_ids,
      pair_recurrence_floor_models = 1 / length(final_ids),
      pair_recurrence_expected_ceiling_bootstrap_uncorrected = 1 - exp(-1)
    ),
    quality_gate_left_fewer_than_two_models = gate$fewer_than_two,
    balance_spread = list(with_gate = spread_with_gate, without_gate = spread_without_gate),
    bootstrap_rank_instability = instability,
    settings = list(source = design$source, formula = formula,
                    matchit_call = design$matchit_call, reference_label =
                      if (is.null(design$reference)) NA_character_ else design$reference$label,
                    algorithms = algorithms, seeds = seeds, folds = folds,
                    cross_fit = cross_fit, n_trees = model_set$settings$n_trees,
                    n_trees_source = model_set$settings$n_trees_source,
                    quality_tolerance = quality_tolerance,
                    balance_rule = balance_rule, balance_margin = balance_margin,
                    balance_cap = balance_cap, balance_threshold = balance_threshold,
                    balance_limit = gate$balance_limit,
                    bootstrap_samples = as.integer(bootstrap_samples),
                    method = "nearest", estimand = "ATT", ratio = 1L,
                    replace = design$replace, m_order = design$m_order,
                    caliper_sd = design$caliper,
                    covariate_names = names(covariates), fold_assignment = model_set$folds,
                    cores = cores, seed = seed),
    components = list(model_set = model_set, match_set = match_set)
  ), class = "propensity_sensitivity")
}

#' @rdname propensity_sensitivity
#' @param digits Significant digits for printed numbers.
#' @param ... Not used.
#' @export
print.propensity_sensitivity <- function(x, digits = 3, ...) {
  value <- function(name, column = "with_gate") {
    table <- x$summary$table
    format_number(table[[column]][table$quantity == name], digits)
  }
  models <- x$models[!x$models$model_id %in% "reference", , drop = FALSE]
  used <- x$summary$models_used
  cat("Propensity model sensitivity\n")
  cat(headline_lines(models, used, x$settings, x$quality_gate_left_fewer_than_two_models,
                     x$summary$table, digits), sep = "\n")
  cat("Mean scaled control spread over treated (descriptive): ",
      value("mean_control_spread_scaled_treated"),
      " (1 = controls differ as much as a typical match differs from its treated patient)\n",
      sep = "")
  used_models <- models[models$model_id %in% used, , drop = FALSE]
  cat("Treated kept: ", format_number(100 * min(used_models$treated_kept), digits),
      " % (lowest model in use)\n", sep = "")
  if (!is.na(x$settings$reference_label)) {
    reference <- x$models[x$models$model_id == "reference", ]
    cat("Reference model (", x$settings$reference_label, "): ",
        if (isTRUE(reference$passed_quality_gate)) "passes the balance gate" else
          "outside the equally-reasonable set",
        "; mean absolute rank difference from the models in use: ",
        value("reference_mean_absolute_rank_difference_treated"), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname propensity_sensitivity
#' @param object A `propensity_sensitivity` object.
#' @export
summary.propensity_sensitivity <- function(object, ...) {
  columns <- c("model_id", "algorithm", "logloss", "smd_max_main", "smd_max_extended",
               "energy_distance", "ks_max", "treated_kept", "passed_quality_gate", "in_use",
               "gate_note")
  structure(list(table = object$summary$table, models = object$models[, columns],
                 quality_gate_left_fewer_than_two_models =
                   object$quality_gate_left_fewer_than_two_models,
                 settings = object$settings),
            class = "summary.propensity_sensitivity")
}

#' @rdname propensity_sensitivity
#' @export
print.summary.propensity_sensitivity <- function(x, digits = 3, ...) {
  models_fitted <- x$models[!x$models$model_id %in% "reference", , drop = FALSE]
  cat(headline_lines(models_fitted, models_fitted$model_id[models_fitted$in_use], x$settings,
                     x$quality_gate_left_fewer_than_two_models, x$table, digits),
      sep = "\n")
  cat("\nStudy-level numbers (with the quality gate / without it):\n")
  table <- x$table
  table$with_gate <- signif(table$with_gate, digits)
  table$without_gate <- signif(table$without_gate, digits)
  print(table, row.names = FALSE)
  cat("\nModels:\n")
  models <- x$models
  for (column in c("logloss", "smd_max_main", "smd_max_extended", "energy_distance", "ks_max",
                   "treated_kept")) {
    models[[column]] <- signif(models[[column]], digits)
  }
  print(models, row.names = FALSE)
  cat("\nQuality gate: log-loss within ", 100 * x$settings$quality_tolerance,
      " % of the best; ", balance_rule_text(x$settings, digits), "\n", sep = "")
  cat("Matching: 1:1 nearest neighbour on the logit, ATT, ",
      if (x$settings$replace) "with" else "without", " replacement, order \"",
      x$settings$m_order, "\", caliper ",
      if (is.null(x$settings$caliper_sd)) "none" else
        paste0(signif(x$settings$caliper_sd, digits), " SD"), "\n", sep = "")
  invisible(x)
}

#' @rdname propensity_sensitivity
#' @param row.names,optional Ignored; present for the generic.
#' @export
# nolint next: object_name_linter. `row.names` is the argument name of the generic.
as.data.frame.propensity_sensitivity <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$patients
}

# ---- internal helpers ------------------------------------------------------------------

unsupported_design_message <- "only 1:1 nearest-neighbour matching is supported in this version"

# Reads the matching design of a matchit object; stops on anything this version cannot redo.
read_matchit_design <- function(x, data, caller) {
  info <- x$info
  if (!identical(info$method, "nearest")) {
    stop(unsupported_design_message, " (the object uses method = \"", info$method, "\").")
  }
  if (is.null(info$ratio) || !isTRUE(all.equal(as.numeric(info$ratio), 1))) {
    stop(unsupported_design_message, " (the object uses ratio = ", info$ratio, ").")
  }
  if (!identical(x$estimand, "ATT")) {
    stop("Only the ATT is supported in this version (the object uses estimand = \"",
         x$estimand, "\").")
  }
  if (is.null(x$distance) || isTRUE(info$mahalanobis)) {
    stop(unsupported_design_message, ", on a propensity score; the object has no ",
         "propensity score distance.")
  }
  if (!is.null(info$discard) && !identical(info$discard, "none")) {
    stop("Objects with `discard` other than \"none\" are not supported in this version.")
  }
  for (argument in c("exact", "mahvars", "antiexact", "s.weights")) {
    if (!is.null(x$call[[argument]])) {
      stop("Objects with `", argument, "` are not supported in this version.")
    }
  }
  n <- length(x$treat)
  if (nrow(data) != n) {
    stop("`data` has ", nrow(data), " rows but the matchit object has ", n, " patients; ",
         "pass the data the object was made with.")
  }
  if (!is.null(names(x$treat)) && !identical(names(x$treat), rownames(data))) {
    stop("The row names of `data` differ from those of the matchit object; pass the data the ",
         "object was made with, in the same row order.")
  }
  caliper <- NULL
  if (length(x$caliper) > 0) {
    caliper_names <- names(x$caliper)
    if (!is.null(caliper_names) && any(nzchar(caliper_names))) {
      stop("Calipers on covariates are not supported in this version; only a caliper on the ",
           "propensity score.")
    }
    caliper <- unname(x$caliper[1]) / stats::sd(x$distance)
  }
  m_order <- tryCatch(eval(x$call$m.order, caller), error = function(e) {
    stop("Could not read `m.order` from the matchit call: ", conditionMessage(e))
  })
  if (is.null(m_order)) m_order <- "largest"
  label <- paste0(info$distance, if (!is.null(info$link)) paste0(", ", info$link))
  list(formula = stats::as.formula(x$formula), replace = isTRUE(info$replace),
       m_order = m_order, caliper = caliper, source = "matchit",
       matchit_call = paste(deparse(x$call, width.cutoff = 500L), collapse = " "),
       reference = list(label = label,
                        rank = rank(as.numeric(x$distance), ties.method = "average") / n,
                        weights = as.numeric(x$weights)))
}

select_covariates <- function(covariates, formula, data) {
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) stop("`covariates` must be a data frame or NULL.")
    if (nrow(covariates) != nrow(data)) {
      stop("`covariates` has ", nrow(covariates), " rows but `data` has ", nrow(data), ".")
    }
    # Non-default row names on both sides must agree: same patients, same order.
    has_row_names <- function(frame) .row_names_info(frame) > 0L
    names_differ <- !identical(rownames(covariates), rownames(data))
    if (has_row_names(covariates) && has_row_names(data) && names_differ) {
      stop("The row names of `covariates` differ from those of `data`; `covariates` must hold ",
           "the same patients in the same order as `data`.")
    }
    return(covariates)
  }
  formula_covariate_columns(formula, data)
}

# The source columns of the right-hand side, resolved with R's own formula rules, so that
# "." is expanded and exclusions such as `- death` are respected. Transformations and
# interactions give their source columns (I(age^2) -> age; age:diabetes -> age, diabetes).
# Only columns of `data` are kept, never the treatment. Order: as in the formula, or as in
# `data` when the formula has ".".
formula_covariate_columns <- function(formula, data) {
  formula_terms <- stats::terms(formula, data = data)
  variables <- unique(unlist(lapply(attr(formula_terms, "term.labels"),
                                    function(label) all.vars(str2lang(label)))))
  right_hand_side <- all.vars(formula[[3]])
  order <- if ("." %in% right_hand_side) names(data) else right_hand_side
  variables <- setdiff(order[order %in% variables], all.vars(formula[[2]]))
  variables <- variables[variables %in% names(data)]
  if (length(variables) == 0) {
    stop("No covariate of the formula was found in `data`.")
  }
  data[, variables, drop = FALSE]
}

# Final quality gate: log-loss part (from the model set) and balance part, then the
# fewer-than-two rule on the final set.
# The balance limit of the final gate: the absolute threshold, or the best refitted
# model's max |SMD| on main terms plus the margin, capped (default) or not ("relative_only").
balance_gate_limit <- function(smd_max_main, rule) {
  if (identical(rule$rule, "absolute")) return(rule$threshold)
  relative_limit <- min(smd_max_main) + rule$margin
  if (identical(rule$rule, "relative_only")) return(relative_limit)
  min(relative_limit, rule$cap)
}

# Shown first by print() and summary() when the gate fell back.
gate_fallback_warning <- paste0("WARNING: fewer than two models passed the quality gate; the ",
                                "diagnostic is outside the setting where it was evaluated.\n",
                                "Every model is used below, so the values do not describe a set ",
                                "of equally good models. Do not read them as a validated ",
                                "result.")

# Names of the algorithms as printed ("bart" -> "BART").
algorithm_display_names <- function(algorithms) {
  labels <- c(logistic_regression = "logistic regression", bart = "BART",
              grf_probability = "GRF", random_forest = "random forest",
              gradient_boosting = "gradient boosting", logistic_splines = "logistic splines",
              lasso = "lasso", lasso_interactions = "lasso with interactions")
  unname(ifelse(algorithms %in% names(labels), labels[algorithms], gsub("_", " ", algorithms)))
}

# The three headline blocks of print() and summary(), in this order: the gate
# composition (models and algorithms kept, why the others were excluded), the study mean of
# rank disagreement, and the balance spread (range of max |SMD| on main terms in use).
headline_lines <- function(models, used, settings, fewer_than_two, table, digits = 3) {
  value <- function(name) format_number(table$with_gate[table$quantity == name], digits)
  in_use <- models[models$model_id %in% used, , drop = FALSE]
  algorithms <- unique(in_use$algorithm)
  failed <- models[!models$passed_quality_gate %in% TRUE, , drop = FALSE]
  reasons <- if (nrow(failed) == 0) character(0) else
    paste0("     ", failed$model_id, ": ", sub("; used anyway.*$", "", failed$gate_note))
  c(if (isTRUE(fewer_than_two)) gate_fallback_warning,
    paste0("1. Quality gate: ", length(used), " of ", nrow(models), " models ",
           if (isTRUE(fewer_than_two)) "in use" else "kept", " (", length(algorithms),
           if (length(algorithms) == 1) " algorithm: " else " algorithms: ",
           paste(algorithm_display_names(algorithms), collapse = ", "), "); log-loss within ",
           100 * settings$quality_tolerance, " % of the best; ",
           balance_rule_text(settings, digits)),
    if (isTRUE(fewer_than_two)) {
      "   Fewer than two models passed, so every model is used. Failed:"
    } else if (nrow(failed) > 0) {
      "   Excluded:"
    },
    reasons,
    paste0("2. Rank disagreement (study mean over treated): ",
           value("mean_rank_disagreement_treated"), "; median between/within ratio ",
           value("median_between_within_ratio_treated"), " (exploratory; near 1 = seed noise)"),
    paste0("   Pair instability across the models in use (study mean over treated; exploratory): ",
           value("pair_instability_models_mean_treated")),
    paste0("   Bootstrap rank instability of the logistic model (study mean over treated): ",
           value("bootstrap_rank_instability_treated")),
    paste0("3. Balance spread (range of max |SMD| main terms across the models in use): ",
           value("smd_max_main_range"), " (", format_number(min(in_use$smd_max_main), digits),
           " to ", format_number(max(in_use$smd_max_main), digits), ")"))
}

# One line stating the balance rule in effect (print and summary).
balance_rule_text <- function(settings, digits = 3) {
  limit <- format_number(settings$balance_limit, digits)
  if (identical(settings$balance_rule, "absolute")) {
    return(paste0("absolute balance rule: max |SMD| main terms <= ", settings$balance_threshold))
  }
  if (identical(settings$balance_rule, "relative_only")) {
    return(paste0("relative-only balance rule (no cap): max |SMD| main terms at most ",
                  settings$balance_margin, " above the best model's, limit here ", limit))
  }
  paste0("relative balance rule: max |SMD| main terms at most ", settings$balance_margin,
         " above the best model's and <= ", settings$balance_cap, ", limit here ", limit)
}

apply_final_quality_gate <- function(models, balance_models, rule) {
  balance <- balance_models[match(models$model_id, balance_models$model_id),
                            c("smd_max_main", "smd_max_extended", "energy_distance", "ks_max")]
  models <- cbind(models, balance)
  balance_limit <- balance_gate_limit(models$smd_max_main, rule)
  balance_reason <- if (identical(rule$rule, "absolute")) {
    "max |SMD| main terms above balance_threshold"
  } else if (identical(rule$rule, "relative_only")) {
    "max |SMD| main terms more than balance_margin above the best model's"
  } else {
    "max |SMD| main terms more than balance_margin above the best model's or above balance_cap"
  }
  models$passed_balance_gate <- models$smd_max_main <= balance_limit
  models$passed_quality_gate <- models$passed_logloss_gate & models$passed_balance_gate
  passed <- models$model_id[models$passed_quality_gate]
  fewer_than_two <- length(passed) < 2
  final_ids <- if (fewer_than_two) models$model_id else passed
  models$in_use <- models$model_id %in% final_ids
  models$is_reference <- FALSE
  models$gate_note <- vapply(seq_len(nrow(models)), function(i) {
    reasons <- c(if (!models$passed_logloss_gate[i])
                   "log-loss above (1 + quality_tolerance) x best",
                 if (!models$passed_balance_gate[i]) balance_reason)
    note <- if (length(reasons) == 0) "passed" else paste(reasons, collapse = "; ")
    if (fewer_than_two && length(reasons) > 0) {
      note <- paste0(note, "; used anyway: fewer than two models passed")
    }
    note
  }, character(1))
  list(models = models, final_ids = final_ids, fewer_than_two = fewer_than_two,
       balance_limit = balance_limit)
}

run_set_diagnostics <- function(model_set, match_set, covariates, ids) {
  list(rank = rank_disagreement(model_set, models = ids),
       control = control_spread(match_set, covariates, models = ids),
       recurrence = pair_recurrence(match_set, models = ids))
}

set_summary <- function(diagnostics, spread, models, ids) {
  rank_summary <- diagnostics$rank$summary
  control_summary <- diagnostics$control$summary
  recurrence_summary <- diagnostics$recurrence$summary
  kept <- models$treated_kept[models$model_id %in% ids]
  balance <- list()
  for (i in seq_len(nrow(spread))) {
    balance[[paste0(spread$statistic[i], "_range")]] <- spread$range[i]
    balance[[paste0(spread$statistic[i], "_sd")]] <- spread$sd[i]
  }
  c(list(mean_rank_disagreement_treated = rank_summary$mean_d_treated,
         mean_rank_disagreement_all = rank_summary$mean_d_all,
         median_between_within_ratio_treated =
           rank_summary$median_between_within_ratio_treated,
         mean_logit_disagreement_treated = rank_summary$mean_d_logit_treated,
         mean_control_spread_scaled_treated = control_summary$mean_scaled_treated,
         p90_control_spread_scaled_treated = control_summary$p90_scaled_treated),
    balance,
    list(pair_recurrence_models_mean_treated = recurrence_summary$mean_models_treated,
         # Pair instability, the second summary of the gated set (exploratory).
         pair_instability_models_mean_treated = 1 - recurrence_summary$mean_models_treated,
         pair_recurrence_bootstrap_mean_treated = recurrence_summary$mean_bootstrap_treated,
         treated_kept_share = mean(kept),
         treated_kept_share_min = min(kept)))
}

reference_model_row <- function(reference, covariates, treatment, balance_limit, columns) {
  statistics <- balance_statistics(covariates, treatment, list(reference = reference$weights))
  passed <- statistics$smd_max_main <= balance_limit
  treated <- treatment == 1
  row <- data.frame(
    model_id = "reference", algorithm = paste0("reference (", reference$label, ")"),
    seed_index = NA_integer_, task_seed = NA_integer_, n_trees = NA_integer_,
    logloss = NA_real_, brier = NA_real_, c_statistic = NA_real_, passed_logloss_gate = NA,
    statistics, passed_balance_gate = passed, passed_quality_gate = passed, in_use = FALSE,
    is_reference = TRUE,
    gate_note = paste0("in-sample score from the user's object, not comparable (log-loss not ",
                       "computed); gate judged on balance only: ",
                       if (passed) "passed" else
                         "reference model outside the equally-reasonable set"),
    treated_kept = sum(reference$weights[treated] > 0) / sum(treated),
    stringsAsFactors = FALSE
  )
  row[, columns]
}

patient_table <- function(diagnostics, model_set, match_set) {
  rank_patients <- diagnostics$rank$patients
  control_patients <- diagnostics$control$patients
  recurrence_patients <- diagnostics$recurrence$patients
  treated <- model_set$treatment == 1
  logistic_score <- rep(NA_real_, length(treated))
  outside_support <- rep(NA, length(treated))
  if ("logistic_regression_1" %in% colnames(model_set$logit)) {
    logistic_score <- unname(model_set$logit[, "logistic_regression_1"])
    control_range <- range(logistic_score[!treated])
    outside_support[treated] <- logistic_score[treated] < control_range[1] |
      logistic_score[treated] > control_range[2]
  }
  data.frame(
    id = match_set$row_ids, treated = as.integer(model_set$treatment),
    rank_disagreement = rank_patients$rank_disagreement,
    rank_disagreement_between = rank_patients$rank_disagreement_between,
    rank_disagreement_within = rank_patients$rank_disagreement_within,
    logit_disagreement = rank_patients$logit_disagreement,
    control_spread = control_patients$control_spread,
    control_spread_scaled = control_patients$control_spread_scaled,
    pair_recurrence_models = recurrence_patients$pair_recurrence_models,
    pair_recurrence_bootstrap = recurrence_patients$pair_recurrence_bootstrap,
    outside_support = outside_support, logistic_score = logistic_score,
    stringsAsFactors = FALSE
  )
}
