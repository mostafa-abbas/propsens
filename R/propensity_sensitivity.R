#' Propensity model sensitivity of a matched design
#'
#' How much does a matched design depend on which of several similarly good propensity models is
#' used? The function refits the propensity model with every algorithm of the model set, repeats
#' the analyst's matching with each model, keeps the models that pass a quality gate (the kept
#' models), and describes how much the kept models disagree. No outcome is used.
#'
#' Every summary is study-level: read it for the study as a whole. Rank disagreement is the
#' primary summary. Pair instability, rank disagreement with equal weight per algorithm and the
#' between/within-algorithm ratio are exploratory. Bootstrap rank instability, balance spread and
#' control spread are descriptive. The per-patient columns of `patients` (and of
#' `as.data.frame()`) describe how the disagreement is spread over the patients; they are not
#' flags for individual patients. No thresholds for a "large" value exist yet.
#'
#' These diagnostics have been studied in the author's simulations and benchmark analyses; a
#' methods paper is in preparation. No peer-reviewed evaluation exists yet.
#'
#' Two situations raise an R warning and a block in `print()` and `summary()`: fewer than two
#' models pass the quality gate (every model is then used, so there is no set of similarly good
#' models to compare), and the kept models are seeds of a single algorithm (rank disagreement and
#' pair instability then measure seed noise within that algorithm).
#'
#' @section Steps:
#' 1. All arguments, the treatment and the covariates are checked before any model is fitted.
#' 2. [fit_propensity_models()] with the formula (cross-fitted scores by default).
#' 3. [match_each_model()] with the matching options and the logistic bootstrap block.
#' 4. [balance_spread()] for every model.
#' 5. Quality gate: `passed_quality_gate = passed_logloss_gate & smd_max_main <=
#'    balance_limit`. With the default `balance_rule = "relative"` the limit is
#'    the smallest max |SMD| on main terms among the refitted models plus `balance_margin`
#'    (0.10), and never above `balance_cap` (0.25): a model is excluded only when its balance is
#'    clearly worse than the best model's. With `balance_rule = "relative_only"` the limit
#'    is the best model's value plus `balance_margin`, without the cap; use it as a sensitivity
#'    analysis when the cap excludes every model (many covariates). With `balance_rule =
#'    "absolute"` the limit is `balance_threshold` (0.10). With `min_treated_kept` set (off by
#'    default), a model that keeps a smaller share of the treated patients fails the gate and
#'    does not enter the best-model value of the relative rule. If fewer than two models pass,
#'    every model is used and `quality_gate_left_fewer_than_two_models` is `TRUE`.
#' 6. [rank_disagreement()], [control_spread()], [pair_recurrence()] and the balance spread on
#'    the kept models (`summary` column `with_gate`) and on every model (`without_gate`).
#'
#' @section Seeds:
#' The fold split and every tree model depend on `seed`. On small data the kept models and every
#' summary can change noticeably with `seed` (on lalonde, rank disagreement varied several-fold
#' between seeds). Rerun with two or three values of `seed` and report all of them;
#' [seed_sensitivity()] does this. Rank disagreement is a standard deviation over the kept
#' models, so it also depends on how many seeds each kept algorithm contributes; "rank
#' disagreement, equal weight per algorithm" counts each algorithm once.
#'
#' @section What the summaries compare under a caliper:
#' With a caliper each model keeps a different subset of the treated patients, shown by
#' `print()` for each kept model. The standardized mean differences compare the matched
#' treated patients with their matched untreated partners; the energy distance
#' (`cobalt::bal.init(estimand = "ATT")`) also measures the distance to the whole treated group,
#' so it also grows when treated patients are dropped. Rank disagreement averages over all
#' treated patients, matched or not; pair instability and control spread use the treated
#' patients matched by at least one kept model. The gate does not look at the share of treated
#' kept unless `min_treated_kept` is set: a model that keeps few treated patients can reach an
#' easy balance and then set the relative limit for the others.
#'
#' @section Estimand:
#' Matching is 1:1 nearest-neighbour on the logit of each model's score, for the average
#' treatment effect in the treated (ATT). Without replacement (the default of the formula
#' interface and the more common applied choice) the result depends on the matching order, and
#' late matches can be poor when untreated patients are not much more numerous than treated
#' ones, so the matched estimate can keep some bias for the ATT. With replacement every treated
#' patient gets its nearest untreated patient; under overlap and a correctly specified
#' propensity model this targets the ATT, at the price of reused untreated patients and a
#' variance that must account for it. With a caliper, treated patients without an eligible
#' partner are dropped and the estimand is the effect in the retained treated patients, not the
#' ATT.
#'
#' @section Reading a matchit object:
#' From a `matchit` object the function reads the formula, the method, `ratio`, `replace`,
#' `m.order` (from the call; MatchIt's default for the ATT with a propensity score is
#' `"largest"`) and the caliper. Only 1:1 nearest-neighbour matching on a propensity score for
#' the ATT is supported, with or without replacement and with or without a caliper on the
#' score; anything else (including `exact`, `mahvars`, `antiexact`, `s.weights`, `reuse.max`,
#' `unit.id` and calipers on covariates) stops with a clear error. MatchIt applies a caliper in
#' standard deviations of the analyst's distance (over all patients); the same number of standard
#' deviations is applied to the logit of every model here. When `data` is `NULL` the data is
#' found by evaluating the object's `data` argument in the calling frame.
#'
#' @section Your own model:
#' The analyst's own score (`x$distance`) is kept as an extra row of `models` with
#' `model_id = "reference"`. Its ranks come from `x$distance` (ranks do not depend
#' on the scale) and its balance statistics from `x$weights`. Its log-loss, Brier score and
#' c-statistic are `NA`: it is an in-sample score from the user's object and is not comparable
#' with the cross-fitted scores. It is judged on balance only: "kept by the gate" or "not kept by
#' the gate (its balance is above the limit)". It is never one of the models compared.
#' `summary` reports its mean absolute rank difference from the kept models, over treated
#' patients: the mean over treated `i` and models `m` of `|r_ref(i) - r_m(i)|`.
#'
#' @param x A `matchit` object (from `MatchIt::matchit()`) or a two-sided formula such as
#'   `treat ~ x1 + x2`. Only baseline covariates belong in the formula; exclude outcomes,
#'   identifiers and post-treatment variables (for example, `treat ~ . - death - id`). The
#'   treatment may be numeric 0/1, logical or a two-level factor (see
#'   [fit_propensity_models()]).
#' @param data The data frame. Required with a formula; with a `matchit` object it is looked up
#'   from the object's call when `NULL`.
#' @param algorithms The propensity model set. The default is logistic regression, BART and the
#'   GRF probability forest (with the default `seeds`, three seeds for each tree algorithm: seven
#'   models; packages 'dbarts' and 'grf'). A faster set is available with
#'   `algorithms = c("logistic_regression", "random_forest", "gradient_boosting")` (packages
#'   'randomForest' and 'gbm'); it has been studied less. See [fit_propensity_models()] for
#'   every option and the fixed settings.
#' @param seeds `1:k`, the number of seeds per tree algorithm (default `1:3`). Only the number
#'   is used; set the randomness with `seed`.
#' @param folds,cross_fit,n_trees,quality_tolerance Passed to
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
#'   instability and the bootstrap pair recurrence: 0 for none, otherwise at least 2
#'   (default 20).
#' @param covariates Data frame of covariates for the control spread and the balance
#'   statistics. Its rows must be the same patients in the same order as `data`; when both
#'   have row names other than the default and they differ, the function stops. `NULL`
#'   (default) uses the source columns of the formula's right-hand side, from `data`, resolved
#'   with R's formula rules: `.` means every other column and exclusions such as `- death`
#'   are respected; `I(age^2)` gives `age`, `age:diabetes` gives `age` and `diabetes`.
#'   Character columns are an error (convert them to factors).
#' @param replace,m_order,caliper Matching options for the formula interface (without replacement,
#'   order `"largest"`, no caliper by default; the caliper is in standard deviations of each
#'   model's logit). With a `matchit` object they are read from the object and must not be
#'   given.
#' @param cores Number of cores (default 1). Run time grows with the number of patients and
#'   covariates (about a minute for 614 patients and 7 covariates; about 16 minutes for 5,735
#'   patients and 50 covariates on one core); on Unix, `cores` runs the models in parallel.
#' @param seed Integer seed from which every other seed is derived (fold split, models,
#'   bootstrap). See the section "Seeds".
#' @param min_treated_kept `NULL` (default, off) or a share between 0 and 1: with a caliper, a
#'   model that keeps a smaller share of the treated patients fails the quality gate and is left
#'   out before the relative balance limit is computed.
#' @param clip_bootstrap `FALSE` (default) or `TRUE`: clip the bootstrap refits' scores to
#'   `[0.001, 0.999]` as the cross-fitted scores are. The default keeps the method as it was
#'   studied; see [match_each_model()].
#' @return An object of class `"propensity_sensitivity"`: a list with
#' * `patients`: one row per patient (`id`, `treated`, `rank_disagreement`,
#'   `rank_disagreement_between`, `rank_disagreement_within`, `logit_disagreement`,
#'   `control_spread`, `control_spread_scaled`, `pair_recurrence_models`,
#'   `pair_recurrence_bootstrap`, `outside_support`, `logistic_score`; the treated-only measures
#'   are `NA` for untreated patients). These columns describe the spread over patients; they
#'   are not per-patient flags.
#' * `models`: one row per model (fit, balance, `treated_kept`, the gate columns, `in_use` =
#'   kept, `gate_note`), plus the row of your own model for a `matchit` object.
#' * `summary`: a list with `table` (the study-level numbers with the gate, `with_gate`, and
#'   without it, `without_gate`), `quality_gate_left_fewer_than_two_models`,
#'   `reference_model_not_kept`, `models_used` (the kept models), `kept_algorithms`,
#'   `kept_algorithms_count`, `kept_models_one_algorithm`, `pair_recurrence_floor_models`,
#'   `pair_instability_maximum` and
#'   `pair_recurrence_expected_ceiling_bootstrap_uncorrected`.
#' * `quality_gate_left_fewer_than_two_models` and `kept_algorithms_count` (also at the top
#'   level, for scripts).
#' * `balance_spread`: spread tables with and without the gate.
#' * `bootstrap_rank_instability`: see [bootstrap_rank_instability()]; `NULL` without
#'   bootstrap refits.
#' * `settings` and `components` (`model_set`, the fitted [fit_propensity_models()] object
#'   with every model's cross-fitted `probability`, `logit` and `rank`, so that the scores can be
#'   reused without refitting; and `match_set`).
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
#' # The default model set with one seed per tree algorithm (three models), to keep the example
#' # short; an analysis uses the default seeds = 1:3 (seven models, about a minute here).
#' propensity_sensitivity(matched, seeds = 1, bootstrap_samples = 5)
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
                                   m_order = "largest", caliper = NULL, cores = 1, seed = 1,
                                   min_treated_kept = NULL, clip_bootstrap = FALSE) {
  # ---- every argument is checked before any model is fitted ----
  algorithms <- match.arg(algorithms, propensity_algorithms, several.ok = TRUE)
  balance_rule <- match.arg(balance_rule)
  for (name in c("balance_margin", "balance_cap", "balance_threshold")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1 || is.na(value) || value < 0) {
      stop("`", name, "` must be a single number of at least 0.")
    }
  }
  check_seeds(seeds)
  check_single_integer(folds, "folds", minimum = 2)
  check_single_integer(cores, "cores", minimum = 1)
  check_single_integer(seed, "seed")
  check_single_integer(bootstrap_samples, "bootstrap_samples", minimum = 0)
  if (bootstrap_samples == 1) {
    stop("`bootstrap_samples` must be 0 (no bootstrap refits) or at least 2: one refit gives ",
         "no spread.")
  }
  check_proportion_or_null(min_treated_kept, "min_treated_kept")
  check_true_false(clip_bootstrap, "clip_bootstrap")
  check_optional_algorithms(algorithms)
  rule <- list(rule = balance_rule, margin = balance_margin, cap = balance_cap,
               threshold = balance_threshold, min_treated_kept = min_treated_kept)
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
  check_matching_options(design$replace, design$m_order, design$caliper)
  formula <- design$formula
  data <- treatment_as_0_1(formula, data)
  prepare_model_data(formula, data)  # treatment, covariates and missing values, checked now
  covariates <- select_covariates(covariates, formula, data)
  check_covariate_frame(covariates)

  # ---- fit, match, balance ----
  model_set <- fit_propensity_models(formula, data, algorithms = algorithms, seeds = seeds,
                                     folds = folds, cross_fit = cross_fit, n_trees = n_trees,
                                     quality_tolerance = quality_tolerance, cores = cores,
                                     seed = seed)
  match_set <- match_each_model(model_set, data, replace = design$replace, ratio = 1,
                                m_order = design$m_order, caliper = design$caliper,
                                bootstrap_samples = bootstrap_samples, cores = cores,
                                seed = seed, clip_bootstrap = clip_bootstrap)
  treatment <- model_set$treatment
  instability <- if (bootstrap_samples >= 2) bootstrap_rank_instability(match_set)
  instability_value <- if (is.null(instability)) NA_real_ else instability$summary$mean_treated
  balance_all <- balance_spread(match_set, treatment, covariates, use_gate = FALSE)

  # ---- quality gate and diagnostics ----
  gate <- apply_final_quality_gate(model_set$models, balance_all$models, rule,
                                   treated_kept = match_set$treated_kept)
  models <- gate$models
  models$treated_kept <- unname(match_set$treated_kept[models$model_id])
  all_ids <- models$model_id
  final_ids <- gate$final_ids
  kept_algorithms <- unique(models$algorithm[models$model_id %in% final_ids])
  one_algorithm <- length(final_ids) >= 2 && length(kept_algorithms) == 1

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

  reference_not_kept <- NA
  if (!is.null(design$reference)) {
    reference_not_kept <- !isTRUE(models$passed_quality_gate[models$model_id == "reference"])
  }
  if (gate$fewer_than_two) warning(fallback_text, call. = FALSE)
  if (one_algorithm) {
    warning(one_algorithm_text(kept_algorithms, gate$fewer_than_two), call. = FALSE)
  }
  structure(list(
    patients = patient_table(with_gate, model_set, match_set),
    models = models,
    summary = list(
      table = summary_table,
      quality_gate_left_fewer_than_two_models = gate$fewer_than_two,
      reference_model_not_kept = reference_not_kept,
      models_used = final_ids,
      kept_algorithms = kept_algorithms,
      kept_algorithms_count = length(kept_algorithms),
      kept_models_one_algorithm = one_algorithm,
      pair_recurrence_floor_models = 1 / length(final_ids),
      pair_instability_maximum = 1 - 1 / length(final_ids),
      pair_recurrence_expected_ceiling_bootstrap_uncorrected = 1 - exp(-1)
    ),
    quality_gate_left_fewer_than_two_models = gate$fewer_than_two,
    kept_algorithms_count = length(kept_algorithms),
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
                    balance_limit = gate$balance_limit, min_treated_kept = min_treated_kept,
                    bootstrap_samples = as.integer(bootstrap_samples),
                    clip_bootstrap = clip_bootstrap,
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
  cat("Propensity model sensitivity\n")
  cat(report_lines(x, digits), sep = "\n")
  invisible(x)
}

#' @rdname propensity_sensitivity
#' @param object A `propensity_sensitivity` object.
#' @export
summary.propensity_sensitivity <- function(object, ...) {
  columns <- c("model_id", "algorithm", "logloss", "smd_max_main", "smd_max_extended",
               "energy_distance", "ks_max", "treated_kept", "passed_quality_gate", "in_use")
  structure(list(table = object$summary$table, models = object$models[, columns],
                 gate_notes = stats::setNames(object$models$gate_note, object$models$model_id),
                 quality_gate_left_fewer_than_two_models =
                   object$quality_gate_left_fewer_than_two_models,
                 settings = object$settings, result = object),
            class = "summary.propensity_sensitivity")
}

#' @rdname propensity_sensitivity
#' @export
print.summary.propensity_sensitivity <- function(x, digits = 3, ...) {
  cat("Propensity model sensitivity: summary\n")
  cat(report_lines(x$result, digits), sep = "\n")
  fallback <- isTRUE(x$quality_gate_left_fewer_than_two_models)
  cat("\nStudy-level numbers",
      if (fallback) " (every model is used, so one column)" else
        " (kept models / every model)", ":\n", sep = "")
  table <- x$table
  shown <- data.frame(quantity = table$quantity,
                      label = unname(summary_labels[table$quantity]),
                      stringsAsFactors = FALSE)
  format_column <- function(values) {
    ifelse(table$quantity == "number_of_models_used", as.character(round(values)),
           vapply(values, format_number, "", digits = digits))
  }
  if (fallback) {
    shown$value <- format_column(table$with_gate)
  } else {
    shown$kept_models <- format_column(table$with_gate)
    shown$every_model <- format_column(table$without_gate)
  }
  print(shown, row.names = FALSE, right = FALSE)
  cat("\nModels (in_use = kept by the quality gate):\n")
  models <- x$models
  for (column in c("logloss", "smd_max_main", "smd_max_extended", "energy_distance", "ks_max",
                   "treated_kept")) {
    models[[column]] <- signif(models[[column]], digits)
  }
  print(models, row.names = FALSE)
  cat("\nGate notes:\n")
  cat(unlist(lapply(names(x$gate_notes), function(id) {
    wrap_text(paste0(id, ": ", x$gate_notes[[id]]), indent = 2, exdent = 6)
  })), sep = "\n")
  cat("\n")
  cat(wrap_text(paste0("Quality gate: log-loss within ", 100 * x$settings$quality_tolerance,
                       " % of the best; ", balance_rule_text(x$settings, digits),
                       if (!is.null(x$settings$min_treated_kept)) {
                         paste0("; share of treated kept at least ",
                                format_setting(x$settings$min_treated_kept))
                       }, ".")), sep = "\n")
  cat(matching_lines(x$result), sep = "\n")
  invisible(x)
}

#' @rdname propensity_sensitivity
#' @param row.names,optional Ignored; present for the generic.
#' @export
# nolint next: object_name_linter. `row.names` is the argument name of the generic.
as.data.frame.propensity_sensitivity <- function(x, row.names = NULL, optional = FALSE, ...) {
  patients <- x$patients
  attr(patients, "note") <- paste(
    "One row per patient. The columns describe how the disagreement is spread over the",
    "patients; they are study-level descriptions, not flags for individual patients."
  )
  patients
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
  for (argument in c("exact", "mahvars", "antiexact", "s.weights", "reuse.max", "unit.id")) {
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

# The balance limit of the final gate: the absolute threshold, or the best refitted
# model's max |SMD| on main terms plus the margin, capped (default) or not ("relative_only").
balance_gate_limit <- function(smd_max_main, rule) {
  if (identical(rule$rule, "absolute")) return(rule$threshold)
  relative_limit <- min(smd_max_main) + rule$margin
  if (identical(rule$rule, "relative_only")) return(relative_limit)
  min(relative_limit, rule$cap)
}

# The R warning and the printed block when the gate fell back.
fallback_text <- paste("WARNING: fewer than two models passed the quality gate, so there is no",
                       "set of similarly good models to compare. Every model is used below;",
                       "read the values as a description of all models, not of a kept set.")

one_algorithm_text <- function(algorithms, fallback = FALSE) {
  paste0("The ", if (fallback) "models in use" else "kept models",
         " are random seeds of one algorithm (", algorithm_display_names(algorithms),
         "): rank disagreement and pair instability here measure seed noise within that ",
         "algorithm, not disagreement between algorithms.")
}

# Names of the algorithms as printed ("bart" -> "BART").
algorithm_display_names <- function(algorithms) {
  labels <- c(logistic_regression = "logistic regression", bart = "BART",
              grf_probability = "GRF", random_forest = "random forest",
              gradient_boosting = "gradient boosting", logistic_splines = "logistic splines",
              lasso = "lasso", lasso_interactions = "lasso with interactions")
  unname(ifelse(algorithms %in% names(labels), labels[algorithms], gsub("_", " ", algorithms)))
}

# Label of each row of the study-level table (see the definitions vignette).
summary_labels <- c(
  number_of_models_used = "count",
  mean_rank_disagreement_treated = "primary",
  mean_rank_disagreement_equal_weight_treated = "exploratory",
  mean_rank_disagreement_all = "descriptive",
  median_between_within_ratio_treated = "exploratory",
  mean_logit_disagreement_treated = "descriptive",
  mean_control_spread_scaled_treated = "descriptive",
  p90_control_spread_scaled_treated = "descriptive",
  smd_max_main_range = "descriptive", smd_max_main_sd = "descriptive",
  smd_max_extended_range = "descriptive", smd_max_extended_sd = "descriptive",
  energy_distance_range = "descriptive", energy_distance_sd = "descriptive",
  ks_max_range = "descriptive", ks_max_sd = "descriptive",
  pair_recurrence_models_mean_treated = "descriptive",
  pair_instability_models_mean_treated = "exploratory",
  pair_recurrence_bootstrap_mean_treated = "descriptive",
  treated_kept_share = "descriptive", treated_kept_share_min = "descriptive",
  bootstrap_rank_instability_treated = "descriptive",
  reference_mean_absolute_rank_difference_treated = "descriptive"
)

# The printed report of print() and summary(), wrapped to the console width.
report_lines <- function(x, digits = 3) {
  table <- x$summary$table
  value <- function(name) format_number(table$with_gate[table$quantity == name], digits)
  models <- x$models[!x$models$model_id %in% "reference", , drop = FALSE]
  used <- x$summary$models_used
  fallback <- isTRUE(x$quality_gate_left_fewer_than_two_models)
  settings <- x$settings
  in_use <- models[models$model_id %in% used, , drop = FALSE]
  algorithms <- unique(in_use$algorithm)
  one_algorithm <- length(used) >= 2 && length(algorithms) == 1
  failed <- models[!models$passed_quality_gate %in% TRUE, , drop = FALSE]
  reasons <- unlist(lapply(seq_len(nrow(failed)), function(i) {
    wrap_text(paste0(failed$model_id[i], ": ",
                     sub("; used anyway.*$", "", failed$gate_note[i])), indent = 5, exdent = 7)
  }))
  kept_shares <- paste0(in_use$model_id, " ",
                        vapply(100 * in_use$treated_kept, format_number, "", digits = digits),
                        " %", collapse = ", ")
  n_used <- length(used)
  algorithm_text <- paste0(length(algorithms),
                           if (length(algorithms) == 1) " algorithm: " else " algorithms: ",
                           paste(algorithm_display_names(algorithms), collapse = ", "))
  share_rule <- if (!is.null(settings$min_treated_kept)) {
    paste0("; share of treated kept at least ", format_setting(settings$min_treated_kept))
  }
  gate_line <- paste0("1. Quality gate: ", n_used, " of ", nrow(models), " models ",
                      if (fallback) "in use" else "kept", " (", algorithm_text,
                      "); log-loss within ", 100 * settings$quality_tolerance, " % of the best; ",
                      balance_rule_text(settings, digits), share_rule, ".")
  c(if (fallback) wrap_text(fallback_text, exdent = 0),
    wrap_text(gate_line),
    if (one_algorithm) wrap_text(one_algorithm_text(algorithms, fallback), indent = 3),
    if (fallback) {
      "   Fewer than two models passed, so every model is used. Failed:"
    } else if (nrow(failed) > 0) {
      "   Excluded:"
    },
    reasons,
    wrap_text(paste0("Treated patients kept by each ", if (fallback) "model" else "kept model",
                     ": ", kept_shares), indent = 3, exdent = 5),
    wrap_text(paste0("2. Rank disagreement (primary; study mean over treated): ",
                     value("mean_rank_disagreement_treated"))),
    wrap_text(paste0("Equal weight per algorithm (exploratory): ",
                     value("mean_rank_disagreement_equal_weight_treated")),
              indent = 3, exdent = 5),
    wrap_text(paste0("Between/within-algorithm ratio (exploratory; median over treated): ",
                     value("median_between_within_ratio_treated"),
                     ". Disagreement between algorithms divided by disagreement between seeds ",
                     "of the same algorithm; with three seeds per algorithm, seed noise alone ",
                     "gives about 0.4."), indent = 3, exdent = 5),
    wrap_text(paste0("Pair instability (exploratory; study mean over treated): ",
                     value("pair_instability_models_mean_treated"), " of a possible ",
                     format_number(1 - 1 / n_used, digits), " with ", n_used, " models. ",
                     "It is high whenever many untreated patients have similar scores, even ",
                     "when the models agree."), indent = 3, exdent = 5),
    wrap_text(paste0("Bootstrap rank instability of the logistic model (descriptive): ",
                     value("bootstrap_rank_instability_treated")), indent = 3, exdent = 5),
    wrap_text(paste0("3. Balance spread (descriptive; range of max |SMD| main terms over the ",
                     if (fallback) "models in use" else "kept models", "): ",
                     value("smd_max_main_range"), " (",
                     format_number(min(in_use$smd_max_main), digits), " to ",
                     format_number(max(in_use$smd_max_main), digits), ")")),
    wrap_text(paste0("Control spread (descriptive; mean scaled over treated): ",
                     value("mean_control_spread_scaled_treated"),
                     "; 1 = the partners chosen by different models differ as much as a ",
                     "typical matched pair (scale: the matches of ", scale_model_name(x), ")")),
    reference_lines(x, digits),
    matching_lines(x))
}

scale_model_name <- function(x) {
  ids <- x$components$match_set$models
  logistic <- ids$model_id[ids$algorithm == "logistic_regression"]
  if (length(logistic) > 0) logistic[1] else x$summary$models_used[1]
}

reference_lines <- function(x, digits = 3) {
  if (is.na(x$settings$reference_label)) return(character(0))
  table <- x$summary$table
  reference <- x$models[x$models$model_id == "reference", ]
  difference <- format_number(table$with_gate[table$quantity ==
                                                "reference_mean_absolute_rank_difference_treated"],
                              digits)
  wrap_text(paste0("Your own model (", x$settings$reference_label, "): ",
                   if (isTRUE(reference$passed_quality_gate)) {
                     "kept by the gate (its balance is within the limit)"
                   } else {
                     "not kept by the gate (its balance is above the limit)"
                   },
                   "; mean absolute rank difference from the ",
                   if (isTRUE(x$quality_gate_left_fewer_than_two_models)) "models in use" else
                     "kept models", ": ", difference))
}

matching_lines <- function(x) {
  settings <- x$settings
  models <- x$models[!x$models$model_id %in% "reference", , drop = FALSE]
  kept <- models$treated_kept[models$model_id %in% x$summary$models_used]
  all_kept <- isTRUE(all(kept >= 1))
  target <- if (all_kept) {
    "target: ATT"
  } else {
    paste0("target: ATT (", format_number(100 * min(kept), 3), " to ",
           format_number(100 * max(kept), 3), " % of treated retained; the estimand is the ",
           "effect in the retained treated patients)")
  }
  wrap_text(paste0("Matching: 1:1 nearest neighbour on the logit, ",
                   if (settings$replace) "with" else "without", " replacement, order \"",
                   settings$m_order, "\", caliper ",
                   if (is.null(settings$caliper_sd)) "none" else
                     paste0(signif(settings$caliper_sd, 3), " SD"), "; ", target, "."))
}

# One line stating the balance rule in effect (print and summary).
balance_rule_text <- function(settings, digits = 3) {
  limit <- format_number(settings$balance_limit, digits)
  if (identical(settings$balance_rule, "absolute")) {
    return(paste0("absolute balance rule: max |SMD| main terms at most ",
                  format_setting(settings$balance_threshold)))
  }
  if (identical(settings$balance_rule, "relative_only")) {
    return(paste0("relative-only balance rule (no cap): max |SMD| main terms at most ",
                  format_setting(settings$balance_margin), " above the best model's (limit here ",
                  limit, ")"))
  }
  paste0("relative balance rule: max |SMD| main terms at most ",
         format_setting(settings$balance_margin), " above the best model's and at most ",
         format_setting(settings$balance_cap), " (limit here ", limit, ")")
}

# Final quality gate: log-loss part (from the model set), the optional share of treated kept,
# and the balance part, then the fewer-than-two rule on the final set.
apply_final_quality_gate <- function(models, balance_models, rule, treated_kept = NULL) {
  balance <- balance_models[match(models$model_id, balance_models$model_id),
                            c("smd_max_main", "smd_max_extended", "energy_distance", "ks_max")]
  models <- cbind(models, balance)
  eligible <- rep(TRUE, nrow(models))
  if (!is.null(rule$min_treated_kept)) {
    eligible <- unname(treated_kept[models$model_id]) >= rule$min_treated_kept
    models$passed_treated_kept <- eligible
  }
  # The best model of the relative rule is taken over the eligible models (all by default).
  balance_limit <- balance_gate_limit(models$smd_max_main[if (any(eligible)) eligible else TRUE],
                                      rule)
  balance_reason <- if (identical(rule$rule, "absolute")) {
    "max |SMD| main terms above balance_threshold"
  } else if (identical(rule$rule, "relative_only")) {
    "max |SMD| main terms more than balance_margin above the best model's"
  } else {
    "max |SMD| main terms more than balance_margin above the best model's or above balance_cap"
  }
  models$passed_balance_gate <- models$smd_max_main <= balance_limit
  models$passed_quality_gate <- models$passed_logloss_gate & models$passed_balance_gate &
    eligible
  passed <- models$model_id[models$passed_quality_gate]
  fewer_than_two <- length(passed) < 2
  final_ids <- if (fewer_than_two) models$model_id else passed
  models$in_use <- models$model_id %in% final_ids
  models$is_reference <- FALSE
  models$gate_note <- vapply(seq_len(nrow(models)), function(i) {
    reasons <- c(if (!models$passed_logloss_gate[i])
                   "log-loss above (1 + quality_tolerance) x best",
                 if (!eligible[i]) "share of treated kept below min_treated_kept",
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
         mean_rank_disagreement_equal_weight_treated = rank_summary$mean_d_between_treated,
         mean_rank_disagreement_all = rank_summary$mean_d_all,
         median_between_within_ratio_treated =
           rank_summary$median_between_within_ratio_treated,
         mean_logit_disagreement_treated = rank_summary$mean_d_logit_treated,
         mean_control_spread_scaled_treated = control_summary$mean_scaled_treated,
         p90_control_spread_scaled_treated = control_summary$p90_scaled_treated),
    balance,
    list(pair_recurrence_models_mean_treated = recurrence_summary$mean_models_treated,
         # Pair instability (exploratory).
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
    gate_note = paste0("your own model: an in-sample score from the matchit object, so its ",
                       "log-loss is not comparable and not computed; judged on balance only: ",
                       if (passed) "kept by the gate (its balance is within the limit)" else
                         "not kept by the gate (its balance is above the limit)"),
    treated_kept = sum(reference$weights[treated] > 0) / sum(treated),
    stringsAsFactors = FALSE
  )
  if ("passed_treated_kept" %in% columns) row$passed_treated_kept <- NA
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
