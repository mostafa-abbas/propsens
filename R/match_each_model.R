#' Match the treated patients once per propensity model
#'
#' Runs 1:1 nearest-neighbour matching for the average treatment effect in the treated (ATT)
#' once per model of a [fit_propensity_models()] object, each time calling
#' `MatchIt::matchit()` with that model's logit column as a numeric `distance`
#' Every matching option is passed explicitly, so a future
#' change of a MatchIt default cannot alter the results. All models are matched, whether or
#' not they passed the log-loss gate; the gate is applied later by the wrapper
#' `propensity_sensitivity()`.
#'
#' @section Matching options:
#' The call is `MatchIt::matchit(formula, data, method = "nearest", estimand = "ATT",
#' distance = <logit column>, replace = replace, ratio = 1, m.order = m_order,
#' caliper = caliper, std.caliper = TRUE, normalize = FALSE)`. The order `"largest"` is
#' MatchIt's default for the ATT when a numeric distance is supplied: the treated patient with
#' the largest score is matched first. A caliper is given in standard deviations of the logit
#' (`std.caliper = TRUE`); it drops treated patients, so the share of treated kept is reported.
#'
#' @section Weights:
#' `normalize = FALSE` is passed so that MatchIt's weights are the raw matching weights: every
#' matched treated patient has weight 1, every matched untreated patient has weight equal to
#' the number of pairs it belongs to, and unmatched patients have weight 0. The sum of the
#' weights of the untreated patients is therefore the number of matched treated patients, also
#' with replacement. (MatchIt's default `normalize = TRUE` rescales those weights to a mean of
#' 1, which changes no weighted mean or balance statistic but would make that sum equal the
#' number of distinct untreated partners.) Each row of a match table carries
#' `weight = 1 / ratio`, which is 1 in this version.
#'
#' @section Bootstrap block:
#' With `bootstrap_samples = B > 0` the function draws `B` nonparametric row bootstrap samples of
#' the data (one derived seed per sample). A sample with fewer than 10 treated or 10 untreated rows
#' is skipped and listed in `bootstrap$skipped`. For each sample, logistic regression with the
#' model set's formula is fitted in-sample, without cross-fitting: main-effects logistic
#' regression barely overfits, and this block measures sampling instability only (see
#' [bootstrap_rank_instability()]). The refitted model predicts logits for every
#' ORIGINAL patient; their percentile ranks form `bootstrap$rank` (n x B). When a factor level
#' of the original data is absent from the sample, the refit uses the full-data design matrix
#' and scores those patients as the reference level (the rule of the cross-fitted logistic
#' model); the samples concerned are listed in
#' `bootstrap$design_rule_samples`. When every level is present but prediction from the
#' refitted model fails (a rank-deficient fit that did not converge), the refit uses the same
#' design-matrix rule; those samples are listed in `bootstrap$prediction_fallback_samples`.
#' No sample is skipped for this reason. The sample is then
#' matched with the same options, rows are mapped back to original ids, and for every
#' original treated patient present in the sample the untreated partner matched to its first
#' occurrence
#' is stored in `bootstrap$pairs` (n x B, `NA` when absent or unmatched). `bootstrap$present`
#' (n x B, logical) says whether each original patient is in each sample. A skipped sample
#' keeps its column, filled with `NA`, so that every matrix is n x B.
#'
#' @section Seeds and parallel runs:
#' Per-task seeds are derived in the main process from `seed`, passed to each worker and
#' applied with [withr::with_seed()]; `m_order = "random"` and the bootstrap draws are
#' therefore reproducible and identical for `cores = 1` and `cores > 1`. Parallel runs use
#' [parallel::mclapply()] on Unix and fall back to a serial run on Windows.
#'
#' @param model_set A `propensity_model_set` from [fit_propensity_models()].
#' @param data The data frame the model set was fitted on (same rows, same order). The
#'   matches are stored with the row names of `data` as ids; the data itself is not stored.
#' @param replace `TRUE` (matching with replacement, the primary analysis) or `FALSE`
#'   (without replacement, the secondary analysis).
#' @param ratio Number of untreated partners per treated patient. Only `1` is supported in this
#'   version.
#' @param m_order Order in which treated patients are matched; one of MatchIt's
#'   `"largest"`, `"smallest"`, `"data"`, `"random"`, `"closest"`, `"farthest"`.
#' @param caliper `NULL` (default, no caliper) or a positive number: the caliper in standard
#'   deviations of the logit.
#' @param oracle_logit Simulation only: a numeric vector of length n with the true logit, for
#'   example the `true_logit` column of [simulate_example_data()]. When given, one extra
#'   entry `oracle` is matched with the same options.
#' @param bootstrap_samples Number of bootstrap samples `B` for the sampling-instability block
#'   (default 0, none).
#' @param cores Number of cores for the matching tasks (default 1).
#' @param seed Integer seed from which every other seed is derived.
#' @param common_random_order Only for `m_order = "random"`: `FALSE` (default) gives each
#'   model its own random order; `TRUE` gives every model (and the oracle) the same random
#'   order, so that the models differ only in their distances, not in the order in which
#'   treated patients are matched. Used by the matching-order sensitivity of the simulation.
#' @return An object of class `"match_set"`: a list with
#' * `matches`: one data frame per model (`treated_id`, `control_id`, `weight`, `pair_id`),
#'   named by model id;
#' * `weights`: n x M matrix of MatchIt weights, rows named by the row ids;
#' * `treated_kept`: named vector, share of treated patients that were matched;
#' * `calls`: named character vector with the deparsed `matchit()` call of each model; the
#'   call names the data and the model set as the caller wrote them and re-evaluates to the
#'   same matches;
#' * `oracle`: `NULL`, or a list with `matches`, `weights`, `treated_kept` and `call`;
#' * `bootstrap`: `NULL`, or a list with `rank`, `pairs`, `present` (n x B matrices) and
#'   `skipped` (data frame of skipped samples with the reason) and `design_rule_samples` (the
#'   samples whose refit used the design-matrix rule for a level absent from the sample) and
#'   `prediction_fallback_samples` (the samples whose refit used it because prediction failed);
#' * `models`: the model table of the model set (`model_id`, `algorithm`, `seed_index`,
#'   `passed_logloss_gate`), used by the diagnostics to apply the quality gate;
#' * `model_ids`, `treatment`, `row_ids`, `settings` and `call`.
#' @examples
#' data("lalonde", package = "MatchIt")
#' model_set <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1, n_trees = 100,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
#' )
#' with_replacement <- match_each_model(model_set, lalonde, replace = TRUE)
#' with_replacement
#' head(with_replacement$matches$logistic_regression_1)
#' without_replacement <- match_each_model(model_set, lalonde, replace = FALSE,
#'                                         bootstrap_samples = 5)
#' without_replacement
#' dim(without_replacement$bootstrap$rank)
#' @importFrom MatchIt matchit
#' @export
match_each_model <- function(model_set, data, replace = TRUE, ratio = 1, m_order = "largest",
                             caliper = NULL, oracle_logit = NULL, bootstrap_samples = 0,
                             cores = 1, seed = 1, common_random_order = FALSE) {
  call <- match.call()
  if (!inherits(model_set, "propensity_model_set")) {
    stop("`model_set` must be a `propensity_model_set` from fit_propensity_models().")
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  n <- length(model_set$treatment)
  if (nrow(data) != n) {
    stop("`data` has ", nrow(data), " rows but the model set was fitted on ", n, " rows. ",
         "Pass the same data frame, in the same row order.")
  }
  data_treatment <- stats::model.response(stats::model.frame(model_set$formula, data = data,
                                                             na.action = stats::na.pass))
  if (!identical(as.integer(data_treatment), as.integer(model_set$treatment))) {
    stop("The treatment in `data` differs from the treatment the model set was fitted on. ",
         "Pass the same data frame, in the same row order.")
  }
  if (!is.logical(replace) || length(replace) != 1 || is.na(replace)) {
    stop("`replace` must be TRUE or FALSE.")
  }
  if (!is.numeric(ratio) || length(ratio) != 1 || is.na(ratio) || ratio != 1) {
    stop("`ratio` other than 1 is not supported in this version.")
  }
  allowed_orders <- c("largest", "smallest", "data", "random", "closest", "farthest")
  if (!is.character(m_order) || length(m_order) != 1 || !m_order %in% allowed_orders) {
    stop("`m_order` must be one of ", paste0('"', allowed_orders, '"', collapse = ", "), ".")
  }
  if (!is.null(caliper)) {
    if (!is.numeric(caliper) || length(caliper) != 1 || is.na(caliper) || caliper <= 0) {
      stop("`caliper` must be NULL or a single positive number (standard deviations of the ",
           "logit).")
    }
  }
  if (!is.null(oracle_logit)) {
    if (!is.numeric(oracle_logit) || length(oracle_logit) != n || anyNA(oracle_logit)) {
      stop("`oracle_logit` must be a numeric vector of length ", n, " without missing values.")
    }
  }
  check_single_integer(bootstrap_samples, "bootstrap_samples", minimum = 0)
  check_single_integer(cores, "cores", minimum = 1)
  check_single_integer(seed, "seed")
  valid_common_order <- is.logical(common_random_order) && length(common_random_order) == 1 &&
    !is.na(common_random_order)
  if (!valid_common_order) stop("`common_random_order` must be TRUE or FALSE.")
  if (common_random_order && m_order != "random") {
    stop("`common_random_order = TRUE` needs `m_order = \"random\"`.")
  }

  row_ids <- rownames(data)
  formula <- model_set$formula
  treatment <- as.integer(model_set$treatment)
  options <- list(replace = replace, ratio = as.integer(ratio), m_order = m_order,
                  caliper = caliper)
  model_ids <- colnames(model_set$logit)
  # The stored calls name the data and the model set the way the caller wrote them.
  data_text <- deparse_one(call$data, "data")
  model_set_text <- deparse_one(call$model_set, "model_set")

  # One matching task per model, plus the oracle; seeds matter only for m_order = "random".
  tasks <- lapply(seq_along(model_ids), function(m) {
    list(name = model_ids[m], logit = model_set$logit[, m],
         distance_text = paste0(model_set_text, '$logit[, "', model_ids[m], '"]'),
         task_seed = derive_seed(seed, if (common_random_order) 100L else 100L + m))
  })
  if (!is.null(oracle_logit)) {
    oracle_seed_index <- 100L + length(model_ids) + 1L
    tasks[[length(tasks) + 1]] <- list(
      name = "oracle", logit = as.numeric(oracle_logit),
      distance_text = deparse_one(call$oracle_logit, "oracle_logit"),
      task_seed = derive_seed(seed, if (common_random_order) 100L else oracle_seed_index)
    )
  }
  match_task <- function(task) {
    match_one_distance(task, formula, data, data_text, treatment, row_ids, options)
  }
  results <- run_tasks(tasks, match_task, cores)
  names(results) <- vapply(tasks, `[[`, "", "name")

  model_results <- results[model_ids]
  weights <- do.call(cbind, lapply(model_results, `[[`, "weights"))
  dimnames(weights) <- list(row_ids, model_ids)

  oracle <- NULL
  if (!is.null(oracle_logit)) {
    oracle <- results[["oracle"]]
    names(oracle$weights) <- row_ids
  }

  bootstrap <- NULL
  if (bootstrap_samples > 0) {
    bootstrap <- bootstrap_logistic_matches(formula, data, treatment, row_ids, options,
                                            bootstrap_samples, seed, cores)
  }

  structure(list(
    matches = lapply(model_results, `[[`, "matches"),
    weights = weights,
    treated_kept = vapply(model_results, `[[`, numeric(1), "treated_kept"),
    calls = vapply(model_results, `[[`, "", "call"),
    oracle = oracle,
    bootstrap = bootstrap,
    model_ids = model_ids,
    models = model_set$models[, c("model_id", "algorithm", "seed_index", "passed_logloss_gate")],
    treatment = treatment,
    row_ids = row_ids,
    settings = list(replace = replace, ratio = as.integer(ratio), m_order = m_order,
                    caliper = caliper, std_caliper = TRUE, normalize = FALSE,
                    estimand = "ATT", method = "nearest", has_oracle = !is.null(oracle_logit),
                    bootstrap_samples = as.integer(bootstrap_samples),
                    common_random_order = common_random_order,
                    bootstrap_minimum_group_size = 10L, cores = cores, seed = seed,
                    rng = rng_settings(), formula = formula),
    call = call
  ), class = "match_set")
}

#' @rdname match_each_model
#' @param x A `match_set` object.
#' @param ... Not used.
#' @export
print.match_set <- function(x, ...) {
  n <- length(x$treatment)
  n_treated <- sum(x$treatment)
  cat("Match set: ", length(x$model_ids), " model(s)", if (!is.null(x$oracle)) " + oracle",
      ", ", n, " patients, ", n_treated, " treated\n", sep = "")
  cat("Matching: 1:", x$settings$ratio, " nearest neighbour on the logit, ATT, ",
      if (x$settings$replace) "with" else "without", " replacement, order \"",
      x$settings$m_order, "\", caliper ",
      if (is.null(x$settings$caliper)) "none" else paste0(x$settings$caliper, " SD"),
      "\n", sep = "")
  entries <- x$matches
  kept <- x$treated_kept
  if (!is.null(x$oracle)) {
    entries$oracle <- x$oracle$matches
    kept <- c(kept, oracle = x$oracle$treated_kept)
  }
  table <- data.frame(
    model_id = names(entries),
    treated_kept = round(kept, 4),
    distinct_controls = vapply(entries, function(m) length(unique(m$control_id)), integer(1)),
    stringsAsFactors = FALSE
  )
  print(table, row.names = FALSE)
  if (!is.null(x$bootstrap)) {
    cat("Bootstrap: ", ncol(x$bootstrap$rank), " logistic refits (",
        nrow(x$bootstrap$skipped), " sample(s) skipped)\n", sep = "")
  }
  invisible(x)
}

# ---- internal helpers ------------------------------------------------------------------

run_tasks <- function(tasks, task_function, cores) {
  use_parallel <- cores > 1 && .Platform$OS.type == "unix"
  if (!use_parallel) return(lapply(tasks, task_function))
  results <- parallel::mclapply(tasks, task_function, mc.cores = cores, mc.preschedule = FALSE,
                                mc.set.seed = FALSE)
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    first_failure <- results[[which(failed)[1]]]
    stop("A matching task failed: ", conditionMessage(attr(first_failure, "condition")))
  }
  results
}

deparse_one <- function(expression, fallback) {
  if (is.null(expression)) return(fallback)
  paste(deparse(expression, width.cutoff = 500L), collapse = " ")
}

build_matchit_call <- function(formula, data_expression, distance_expression, options) {
  as.call(list(
    quote(MatchIt::matchit), formula = formula, data = data_expression, method = "nearest",
    estimand = "ATT", distance = distance_expression, replace = options$replace,
    ratio = options$ratio, m.order = options$m_order, caliper = options$caliper,
    std.caliper = TRUE, normalize = FALSE
  ))
}

run_matchit <- function(formula, data, distance, options, task_seed) {
  matchit_call <- build_matchit_call(formula, quote(data), quote(distance), options)
  with_task_seed(task_seed, eval(matchit_call))
}

# Turns a matchit object into the stored pieces: match table, weights and share of treated kept.
read_matchit_object <- function(matched, treatment, row_ids, options) {
  match_matrix <- matched$match.matrix
  treated_rows <- match(rownames(match_matrix), row_ids)
  control_rows <- match(match_matrix[, 1], row_ids)
  kept <- !is.na(control_rows)
  matches <- data.frame(
    treated_id = row_ids[treated_rows[kept]],
    control_id = row_ids[control_rows[kept]],
    weight = rep(1 / options$ratio, sum(kept)),
    pair_id = seq_len(sum(kept)),
    stringsAsFactors = FALSE
  )
  weights <- as.numeric(matched$weights)
  names(weights) <- row_ids
  list(matches = matches, weights = weights,
       treated_kept = sum(kept) / sum(treatment == 1),
       treated_rows = treated_rows, control_rows = control_rows)
}

match_one_distance <- function(task, formula, data, data_text, treatment, row_ids, options) {
  matched <- run_matchit(formula, data, task$logit, options, task$task_seed)
  pieces <- read_matchit_object(matched, treatment, row_ids, options)
  stored_call <- build_matchit_call(formula, str2lang(data_text), str2lang(task$distance_text),
                                    options)
  list(matches = pieces$matches, weights = pieces$weights, treated_kept = pieces$treated_kept,
       call = paste(deparse(stored_call, width.cutoff = 500L), collapse = " "))
}

bootstrap_logistic_matches <- function(formula, data, treatment, row_ids, options,
                                       bootstrap_samples, seed, cores) {
  n <- length(treatment)
  minimum_group <- 10L
  tasks <- lapply(seq_len(bootstrap_samples), function(b) {
    list(index = b, draw_seed = derive_seed(seed, 1000L + b),
         match_seed = derive_seed(seed, 2000L + b))
  })
  one_sample <- function(task) {
    rows <- with_task_seed(task$draw_seed, sample.int(n, n, replace = TRUE))
    n_treated <- sum(treatment[rows] == 1)
    n_control <- sum(treatment[rows] == 0)
    if (n_treated < minimum_group || n_control < minimum_group) {
      return(list(skipped = paste0("fewer than ", minimum_group, " treated or controls (",
                                   n_treated, " treated, ", n_control, " controls)")))
    }
    resample <- data[rows, , drop = FALSE]
    rownames(resample) <- NULL
    logits <- bootstrap_sample_logits(formula, data, resample, rows, prepared)
    if (!is.null(logits$skipped)) return(list(skipped = logits$skipped))
    original_logit <- logits$original_logit
    resample_logit <- logits$resample_logit
    matched <- run_matchit(formula, resample, resample_logit, options, task$match_seed)
    pieces <- read_matchit_object(matched, treatment[rows], as.character(seq_len(n)), options)
    # Map resample rows back to original rows; keep the first occurrence of each treated.
    original_treated <- rows[pieces$treated_rows]
    original_control <- rows[pieces$control_rows]
    first_occurrence <- !duplicated(original_treated)
    pairs <- rep(NA_character_, n)
    pairs[original_treated[first_occurrence]] <-
      ifelse(is.na(original_control[first_occurrence]), NA_character_,
             row_ids[original_control[first_occurrence]])
    present <- logical(n)
    present[rows] <- TRUE
    list(rank = rank(original_logit, ties.method = "average") / n, pairs = pairs,
         present = present, path = logits$path)
  }
  prepared <- prepare_model_data(formula, data)
  results <- run_tasks(tasks, one_sample, cores)
  skipped_flag <- vapply(results, function(r) !is.null(r$skipped), logical(1))
  skipped <- data.frame(
    sample = which(skipped_flag),
    reason = vapply(results[skipped_flag], `[[`, "", "skipped"),
    stringsAsFactors = FALSE
  )
  # Skipped samples keep their column (all NA) so the matrices stay n x B.
  bind <- function(field, empty) {
    columns <- lapply(results, function(r) if (is.null(r$skipped)) r[[field]] else rep(empty, n))
    out <- do.call(cbind, columns)
    dimnames(out) <- list(row_ids, paste0("bootstrap_", seq_len(bootstrap_samples)))
    out
  }
  path <- vapply(results, function(r) if (is.null(r$path)) NA_character_ else r$path, "")
  list(rank = bind("rank", NA_real_), pairs = bind("pairs", NA_character_),
       present = bind("present", NA), skipped = skipped,
       design_rule_samples = which(path %in% "design_unseen_levels"),
       prediction_fallback_samples = which(path %in% "design_fallback_predict_error"))
}

# Logits of the original data from a bootstrap refit (a separate function so that a test can
# make it fail as a non-converged rank-deficient fit does on real data).
predict_bootstrap_refit <- function(fit, data) {
  as.numeric(stats::predict(fit, newdata = data, type = "link"))
}

# Logits of one bootstrap refit of the logistic model: for every original patient
# (`original_logit`) and for the resample rows (`resample_logit`), or `skipped` with the reason.
# When a factor level of the original data is absent from the resample,
# predict() on the original data would stop ("factor has new levels"); the refit then uses the
# design-matrix rule of the cross-fitted logistic model (logistic_on_design(): full-data dummy
# columns, a coefficient that cannot be estimated set to 0, so those patients are scored as the
# reference level). When the formula refit succeeds but predict() on the
# original data stops (a rank-deficient fit that did not converge can leave non-finite values in
# its QR), the refit also uses the design-matrix rule, as the cross-fitted logistic model does
# ("design_fallback_predict_error"). Otherwise the formula path is unchanged. `path` names the
# path used. `unseen_levels_rule = "skip"` and `predict_error_rule = "skip"` give the earlier
# behaviour (the sample is skipped); they are kept only to check saved results.
bootstrap_sample_logits <- function(formula, data, resample, rows, prepared,
                                    unseen_levels_rule = c("design", "skip"),
                                    predict_error_rule = c("design", "skip")) {
  unseen_levels_rule <- match.arg(unseen_levels_rule)
  predict_error_rule <- match.arg(predict_error_rule)
  on_design <- function(path) {
    design <- cbind(1, prepared$design)
    fit <- tryCatch(stats::glm.fit(design[rows, , drop = FALSE], prepared$treatment[rows],
                                   family = stats::binomial()),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      return(list(skipped = paste0("logistic refit failed: ", conditionMessage(fit))))
    }
    coefficients <- fit$coefficients
    coefficients[is.na(coefficients)] <- 0
    original_logit <- as.numeric(design %*% coefficients)
    list(original_logit = original_logit, resample_logit = original_logit[rows],
         design_rule = path == "design_unseen_levels", path = path)
  }
  if (unseen_levels_rule == "design" &&
        has_unseen_levels(formula, data, rows, seq_len(nrow(data)))) {
    return(on_design("design_unseen_levels"))
  }
  fit <- tryCatch(stats::glm(formula, data = resample, family = stats::binomial()),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(skipped = paste0("logistic refit failed: ", conditionMessage(fit))))
  }
  original_logit <- tryCatch(predict_bootstrap_refit(fit, data), error = function(e) e)
  if (inherits(original_logit, "error")) {
    # A level absent from the resample under the "skip" rule also stops here; only a failure
    # with every level present is a prediction failure for the design-matrix fallback.
    if (predict_error_rule == "design" &&
          !has_unseen_levels(formula, data, rows, seq_len(nrow(data)))) {
      return(on_design("design_fallback_predict_error"))
    }
    return(list(skipped = paste0("prediction failed: ", conditionMessage(original_logit))))
  }
  list(original_logit = original_logit,
       resample_logit = as.numeric(stats::predict(fit, type = "link")), design_rule = FALSE,
       path = "formula")
}
