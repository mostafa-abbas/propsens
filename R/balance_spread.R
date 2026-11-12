#' Study-level balance spread across propensity models
#'
#' Descriptive. For each model, the model's MatchIt weights from
#' [match_each_model()] are passed to cobalt and four balance statistics are computed:
#' * `smd_max_main`: maximum absolute standardized mean difference over the main terms;
#' * `smd_max_extended`: the same over main terms plus squares and 2-way interactions
#'   (`cobalt::bal.tab(int = interactions, poly = polynomial)`);
#' * `energy_distance`: `cobalt::bal.compute()` with `stat = "energy.dist"`;
#' * `ks_max`: `cobalt::bal.compute()` with `stat = "ks.max"`.
#'
#' With many covariates `smd_max_extended` is usually driven by sparse squares and interactions
#' (a rare level times another covariate) and can be large even when every main term is well
#' balanced; report the main-terms range (`smd_max_main`) as the headline balance spread. With a
#' caliper the standardized mean differences compare the matched treated patients with their
#' matched untreated partners, while the energy distance (`estimand = "ATT"`) also measures the
#' distance to the whole treated group, so it also grows when treated patients are dropped.
#'
#' The spread of each statistic over the kept models is its range (max - min), its standard
#' deviation, and the worst (largest) value with the model that produced it. The table of
#' models lists every model of the match set with an `in_use` column, so that the balance part
#' of the quality gate can read every model; the spread uses the kept models only. The
#' oracle entry and the bootstrap block are never used.
#'
#' @section cobalt settings (verified against cobalt 5.0.0):
#' `estimand = "ATT"`, `method = "matching"`, `s.d.denom = "treated"`: the standardization uses
#' the standard deviation of the treated group in the unweighted sample. `binary = "std"` and
#' `continuous = "std"`: binary covariates are standardized as well (cobalt's `bal.tab()`
#' default for binary covariates is the raw difference in proportions, while
#' `bal.compute(stat = "smd.max")` standardizes them; the package uses the standardized version
#' throughout so the two agree). Example: on lalonde after the logistic
#' match with replacement, `married` has a raw difference in proportions of 0.059 and a
#' standardized difference of 0.152, so the choice decides whether a 0.10 balance threshold is
#' met. Weighted means, the energy distance and the KS statistics do
#' not depend on how the weights of the untreated patients are scaled.
#'
#' @inheritSection rank_disagreement Quality gate
#' @param match_set A `match_set` from [match_each_model()].
#' @param treatment The 0/1 treatment vector, in the data's row order; must equal the treatment
#'   stored in `match_set`.
#' @param covariates A data frame with one row per patient, in the data's row order. Factors are
#'   split into dummies by cobalt; character columns are an error.
#' @param interactions `TRUE` (default) to include all 2-way interactions in `smd_max_extended`.
#' @param polynomial Highest power of the continuous covariates in `smd_max_extended`
#'   (default 2, squares).
#' @param use_gate `TRUE` (default) to compute the spread over the models that passed the
#'   log-loss gate; `FALSE` to use every model.
#' @param models Optional character vector of model ids to use instead of a gate; it
#'   overrides `use_gate`. The wrapper [propensity_sensitivity()] passes the final
#'   quality-gate set (log-loss and balance) this way.
#' @return An object of class `"balance_spread"`: a list with `models` (one row per model:
#'   `model_id`, `in_use`, `smd_max_main`, `smd_max_extended`, `energy_distance`, `ks_max`),
#'   `spread` (one row per statistic: `statistic`, `range`, `sd`, `worst_value`, `worst_model`),
#'   `summary`, `models_used`, `quality_gate_left_fewer_than_two_models` and `settings`.
#' @examples
#' data("lalonde", package = "MatchIt")
#' model_set <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1, n_trees = 100,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
#' )
#' match_set <- match_each_model(model_set, lalonde, replace = TRUE)
#' covariates <- lalonde[, c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
#' balance <- balance_spread(match_set, lalonde$treat, covariates)
#' balance
#' @importFrom cobalt bal.tab bal.init bal.compute
#' @export
balance_spread <- function(match_set, treatment, covariates, interactions = TRUE,
                           polynomial = 2, use_gate = TRUE, models = NULL) {
  check_match_set(match_set)
  check_true_false(interactions, "interactions")
  check_true_false(use_gate, "use_gate")
  check_single_integer(polynomial, "polynomial", minimum = 1)
  n <- length(match_set$treatment)
  if (!is.numeric(treatment) || length(treatment) != n ||
        !identical(as.integer(treatment), as.integer(match_set$treatment))) {
    stop("`treatment` must be the 0/1 treatment the match set was built on (length ", n,
         ", same row order).")
  }
  if (!is.data.frame(covariates)) stop("`covariates` must be a data frame.")
  if (nrow(covariates) != n) {
    stop("`covariates` has ", nrow(covariates), " rows but the match set has ", n,
         " patients. Pass the covariates in the data's row order.")
  }
  character_columns <- names(covariates)[vapply(covariates, is.character, logical(1))]
  if (length(character_columns) > 0) {
    stop("Character column(s) in `covariates`: ", paste(character_columns, collapse = ", "),
         ". Convert them to factors first.")
  }
  if (anyNA(covariates)) stop("Missing values in `covariates`. Impute first.")

  in_use <- models_in_use(match_set, use_gate, models)
  all_ids <- match_set$model_ids
  treatment <- as.integer(treatment)
  weight_list <- lapply(all_ids, function(id) unname(match_set$weights[, id]))
  names(weight_list) <- all_ids

  statistics_table <- balance_statistics(covariates, treatment, weight_list, interactions,
                                         polynomial)
  models <- data.frame(model_id = all_ids, in_use = all_ids %in% in_use$model_ids,
                       statistics_table, stringsAsFactors = FALSE)
  spread <- balance_spread_table(models, in_use$model_ids)
  summary <- stats::setNames(as.list(spread$range), paste0(spread$statistic, "_range"))
  structure(list(
    models = models, spread = spread, summary = summary,
    models_used = in_use$model_ids,
    quality_gate_left_fewer_than_two_models = in_use$quality_gate_left_fewer_than_two_models,
    settings = list(use_gate = use_gate, gate_applied = in_use$gate_applied,
                    models_supplied = in_use$models_supplied,
                    interactions = interactions, polynomial = polynomial, estimand = "ATT",
                    s.d.denom = "treated", binary = "std", continuous = "std")
  ), class = "balance_spread")
}

#' @rdname balance_spread
#' @param x A `balance_spread` object.
#' @param ... Not used.
#' @export
print.balance_spread <- function(x, ...) {
  print_diagnostic_header("Balance spread", x)
  table <- x$models
  for (column in c("smd_max_main", "smd_max_extended", "energy_distance", "ks_max")) {
    table[[column]] <- signif(table[[column]], 4)
  }
  print(table, row.names = FALSE)
  cat("Spread over the models in use:\n")
  spread <- x$spread
  spread$range <- signif(spread$range, 4)
  spread$sd <- signif(spread$sd, 4)
  spread$worst_value <- signif(spread$worst_value, 4)
  print(spread, row.names = FALSE)
  invisible(x)
}

# The four balance statistics for each entry of a named list of weight vectors:
# maximum |SMD| on main terms, maximum |SMD| with powers up to `polynomial` and 2-way
# interactions, energy distance, maximum KS statistic. Settings as in the help page (section
# "cobalt settings").
balance_statistics <- function(covariates, treatment, weight_list, interactions = TRUE,
                               polynomial = 2) {
  ids <- names(weight_list)
  max_abs_smd <- function(int, poly, columns = names(covariates)) {
    table <- cobalt::bal.tab(covariates[, columns, drop = FALSE], treat = treatment,
                             weights = weight_list, method = "matching", estimand = "ATT",
                             s.d.denom = "treated",
                             binary = "std", continuous = "std", int = int, poly = poly,
                             un = FALSE, quick = TRUE)$Balance
    # cobalt names the column "Diff.<name>" for a list of several weights, "Diff.Adj" for one.
    column_of <- function(id) {
      if (length(ids) == 1) "Diff.Adj" else paste0("Diff.", id)
    }
    vapply(ids, function(id) max(abs(table[[column_of(id)]]), na.rm = TRUE), numeric(1))
  }
  compute_stat <- function(stat) {
    init <- cobalt::bal.init(covariates, treat = treatment, stat = stat, estimand = "ATT")
    vapply(ids, function(id) cobalt::bal.compute(init, weights = weight_list[[id]]),
           numeric(1))
  }
  data.frame(smd_max_main = unname(max_abs_smd(FALSE, 1)),
             smd_max_extended = unname(extended_max_abs_smd(covariates, max_abs_smd,
                                                            interactions, polynomial)),
             energy_distance = unname(compute_stat("energy.dist")),
             ks_max = unname(compute_stat("ks.max")))
}

# Maximum |SMD| over main terms, powers and 2-way interactions without building every term
# at once. cobalt builds each power and each product from its own components and drops a term by
# its own values only, so every term of the full table appears, with the same SMD, in the table
# of some pair of covariate blocks; the maximum over block pairs is the full maximum. With at
# most `block_size` covariates there is one block and the full table is used, as before.
extended_max_abs_smd <- function(covariates, max_abs_smd, interactions, polynomial,
                                 block_size = 20L) {
  columns <- names(covariates)
  if (length(columns) <= block_size) return(max_abs_smd(interactions, polynomial))
  blocks <- split(columns, ceiling(seq_along(columns) / block_size))
  pairs <- utils::combn(length(blocks), 2L, simplify = FALSE)
  maxima <- lapply(pairs, function(pair) {
    max_abs_smd(interactions, polynomial, c(blocks[[pair[1]]], blocks[[pair[2]]]))
  })
  Reduce(pmax, maxima)
}

# Spread of the four statistics over the models `used_ids` of a per-model table: range, standard
# deviation (n - 1; NA for one model), worst value and the model that has it.
balance_spread_table <- function(models, used_ids) {
  statistics <- c("smd_max_main", "smd_max_extended", "energy_distance", "ks_max")
  used <- models[models$model_id %in% used_ids, , drop = FALSE]
  do.call(rbind, lapply(statistics, function(statistic) {
    values <- used[[statistic]]
    worst <- which.max(values)
    data.frame(statistic = statistic, range = max(values) - min(values),
               sd = if (length(values) > 1) stats::sd(values) else NA_real_,
               worst_value = values[worst], worst_model = used$model_id[worst],
               stringsAsFactors = FALSE)
  }))
}
