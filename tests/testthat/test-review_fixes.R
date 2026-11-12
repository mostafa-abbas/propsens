# One test per behaviour added in the review of 0.1.0: warnings, early checks, unsupported
# matchit options, treatment coding, bootstrap messages and options, the seed helper.
data("lalonde", package = "MatchIt")
review_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75
review_fast <- c("logistic_regression", "random_forest", "gradient_boosting")
quick_run <- function(...) {
  # nolint next: object_usage_linter. `lalonde` is loaded by data() at the top of this file.
  propensity_sensitivity(review_formula, data = lalonde, algorithms = review_fast, seeds = 1,
                         n_trees = 50, bootstrap_samples = 0, ...)
}

# Captures the warnings of an expression and returns its value with the messages attached.
collect_warnings <- function(expression) {
  messages <- character(0)
  value <- withCallingHandlers(expression, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = messages)
}

test_that("a kept set from one algorithm raises a warning, a flag and a printed block", {
  # Only random forests: every kept model is a seed of one algorithm.
  run <- collect_warnings(propensity_sensitivity(review_formula, data = lalonde,
                                                 algorithms = "random_forest", seeds = 1:2,
                                                 bootstrap_samples = 0, replace = TRUE))
  result <- run$value
  expect_true(any(grepl("are random seeds of one algorithm (random forest)", run$warnings,
                        fixed = TRUE)))
  expect_identical(result$kept_algorithms_count, 1L)
  expect_identical(result$summary$kept_algorithms_count, 1L)
  printed <- utils::capture.output(print(result))
  gate <- grep("^1\\. Quality gate:", printed)
  block <- grep("random seeds of one algorithm", printed)
  expect_length(block, 1)
  expect_gt(block, gate)
  expect_match(printed_text(result),
               "measure seed noise within that algorithm, not disagreement between algorithms")
  # With two algorithms kept, no such warning.
  two <- collect_warnings(quick_run(replace = TRUE, seed = 2))
  if (two$value$kept_algorithms_count >= 2) {
    expect_false(any(grepl("random seeds of one algorithm", two$warnings)))
  }
})

test_that("the gate fallback raises an R warning as well as the printed block", {
  run <- collect_warnings(quick_run(balance_rule = "absolute", balance_threshold = 0.001))
  expect_true(run$value$quality_gate_left_fewer_than_two_models)
  fallback_start <- paste("WARNING: fewer than two models passed the quality gate, so there is",
                          "no set of similarly good models to compare.")
  expect_true(any(startsWith(run$warnings, fallback_start)))
  expect_true(any(grepl("^WARNING", utils::capture.output(print(run$value)))))
})

test_that("arguments are checked before any model is fitted", {
  local_mocked_bindings(fit_propensity_models = function(...) stop("a model was fitted"))
  expect_error(quick_run(caliper = -1), "`caliper` must be NULL or a single positive number")
  expect_error(quick_run(m_order = "biggest"), "`m_order` must be one of")
  expect_error(quick_run(replace = NA), "`replace` must be TRUE or FALSE")
  expect_error(propensity_sensitivity(review_formula, data = lalonde, algorithms = review_fast,
                                      seeds = 1, bootstrap_samples = 1),
               "`bootstrap_samples` must be 0 .* or at least 2")
  expect_error(quick_run(min_treated_kept = 2), "`min_treated_kept`")
  expect_error(quick_run(clip_bootstrap = "yes"), "`clip_bootstrap`")
  with_character <- lalonde
  with_character$race <- as.character(with_character$race)
  expect_error(propensity_sensitivity(review_formula, data = with_character,
                                      algorithms = review_fast, seeds = 1),
               "Character column\\(s\\) in `covariates`: race")
  analyst <- MatchIt::matchit(review_formula, data = lalonde, replace = TRUE)
  expect_error(propensity_sensitivity(analyst, data = lalonde, seeds = 2:3), "`seeds` must be 1:k")
})

test_that("matchit objects with reuse.max or unit.id are rejected", {
  with_reuse <- MatchIt::matchit(review_formula, data = lalonde, replace = TRUE, reuse.max = 3)
  expect_error(propensity_sensitivity(with_reuse, data = lalonde),
               "Objects with `reuse.max` are not supported")
  with_unit <- lalonde
  with_unit$unit <- rep(seq_len(nrow(lalonde) / 2), each = 2)
  unit_object <- MatchIt::matchit(review_formula, data = with_unit, unit.id = ~unit)
  expect_error(propensity_sensitivity(unit_object, data = with_unit),
               "Objects with `unit.id` are not supported")
})

test_that("`seeds` must be 1:k; only the number of seeds is used", {
  expect_error(fit_propensity_models(review_formula, lalonde, algorithms = review_fast,
                                     seeds = c(101, 999), n_trees = 50),
               "only the number of seeds is used; set the randomness with `seed`")
  expect_error(propensity_sensitivity(review_formula, data = lalonde, seeds = 0),
               "`seeds` must be 1:k")
  expect_silent(check_seeds(1:3))
  expect_silent(check_seeds(c(1, 2)))
})

test_that("a logical or two-level factor treatment gives the same result as 0/1", {
  numeric_run <- quiet_gate(quick_run(replace = TRUE, seed = 4))
  logical_data <- lalonde
  logical_data$treat <- lalonde$treat == 1
  logical_run <- quiet_gate(propensity_sensitivity(review_formula, data = logical_data,
                                                   algorithms = review_fast, seeds = 1,
                                                   n_trees = 50, bootstrap_samples = 0,
                                                   replace = TRUE, seed = 4))
  factor_data <- lalonde
  factor_data$treat <- factor(ifelse(lalonde$treat == 1, "treated", "control"))
  factor_run <- quiet_gate(propensity_sensitivity(review_formula, data = factor_data,
                                                  algorithms = review_fast, seeds = 1,
                                                  n_trees = 50, bootstrap_samples = 0,
                                                  replace = TRUE, seed = 4))
  expect_equal(logical_run$patients, numeric_run$patients)
  expect_equal(factor_run$patients, numeric_run$patients)
  expect_equal(factor_run$summary$table, numeric_run$summary$table)
  # MatchIt's rule: a factor with levels "0"/"1" or numbers codes the smaller (or 0) as untreated;
  # other levels code the first level as untreated.
  expect_identical(coerce_treatment(factor(c("yes", "no", "yes"))), c(1L, 0L, 1L))
  expect_identical(coerce_treatment(factor(c("2", "5", "5"))), c(0L, 1L, 1L))
  expect_identical(coerce_treatment(c(TRUE, FALSE)), c(1L, 0L))
  # A matchit object made with a factor treatment is read as well.
  analyst <- MatchIt::matchit(review_formula, data = factor_data, replace = TRUE)
  expect_s3_class(quiet_gate(propensity_sensitivity(analyst, data = factor_data,
                                                    algorithms = review_fast, seeds = 1,
                                                    n_trees = 50, bootstrap_samples = 0)),
                  "propensity_sensitivity")
})

test_that("skipped bootstrap samples are reported in a message", {
  set.seed(3)
  small <- data.frame(x1 = stats::rnorm(60), x2 = stats::rnorm(60))
  small$treat <- c(rep(1L, 7), rep(0L, 53))
  small_set <- fit_propensity_models(treat ~ x1 + x2, small, algorithms = "logistic_regression",
                                     cross_fit = FALSE)
  expect_message(match_each_model(small_set, small, replace = TRUE, bootstrap_samples = 4),
                 "Bootstrap: 4 of 4 samples skipped \\(fewer than 10 treated or untreated: 4\\)")
})

test_that("the scaled control spread is NA, not infinite, when the scale distance is 0", {
  model_set <- fit_propensity_models(review_formula, lalonde, algorithms = review_fast,
                                     seeds = 1, n_trees = 50)
  match_set <- match_each_model(model_set, lalonde, replace = FALSE)
  covariates <- lalonde[, c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
  logistic_pairs <- match_set$matches$logistic_regression_1
  # Make every logistic partner identical to its treated patient: the scale distance is 0.
  same <- covariates
  same[match(logistic_pairs$control_id, rownames(lalonde)), ] <-
    covariates[match(logistic_pairs$treated_id, rownames(lalonde)), ]
  spread <- suppressMessages(control_spread(match_set, same, use_gate = FALSE))
  expect_equal(spread$summary$scale_distance, 0)
  expect_false(any(is.infinite(spread$patients$control_spread_scaled)))
  expect_true(all(is.na(spread$patients$control_spread_scaled)))
  expect_true(is.na(spread$summary$mean_scaled_treated))
  expect_output(print(spread), "logistic_regression_1 = 0 \\(scaled values set to NA\\)")
})

test_that("min_treated_kept (off by default) removes models that keep few treated patients", {
  log_loss_table <- data.frame(model_id = c("a", "b", "c"),
                               passed_logloss_gate = c(TRUE, TRUE, TRUE))
  balance_table <- data.frame(model_id = c("a", "b", "c"), smd_max_main = c(0.00, 0.15, 0.20),
                              smd_max_extended = 0, energy_distance = 0, ks_max = 0)
  kept <- c(a = 0.02, b = 0.80, c = 0.85)
  rule <- list(rule = "relative", margin = 0.10, cap = 0.25)
  # Default: model a (2 % kept, perfect balance) sets the limit 0.10, so b and c fail.
  default_gate <- apply_final_quality_gate(log_loss_table, balance_table, rule, kept)
  expect_equal(default_gate$balance_limit, 0.10)
  expect_true(default_gate$fewer_than_two)
  expect_false("passed_treated_kept" %in% names(default_gate$models))
  # With the option: a is left out first; the limit comes from b (0.25) and b, c pass.
  rule$min_treated_kept <- 0.5
  gate <- apply_final_quality_gate(log_loss_table, balance_table, rule, kept)
  expect_equal(gate$balance_limit, 0.25)
  expect_identical(gate$final_ids, c("b", "c"))
  expect_match(gate$models$gate_note[1], "share of treated kept below min_treated_kept")
  # End to end the option is recorded and shown.
  result <- quiet_gate(quick_run(caliper = 0.05, min_treated_kept = 0.9))
  expect_equal(result$settings$min_treated_kept, 0.9)
  fitted <- result$models
  expect_false(any(fitted$passed_quality_gate[fitted$treated_kept < 0.9]))
  expect_match(printed_text(result), "share of treated kept at least 0.90")
})

test_that("clip_bootstrap clips the bootstrap scores (off by default)", {
  set.seed(5)
  n <- 300
  separated <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  separated$treat <- as.integer(separated$x1 > 0.5)
  set <- suppressWarnings(fit_propensity_models(treat ~ x1 + x2, separated,
                                                algorithms = "logistic_regression"))
  plain <- suppressWarnings(match_each_model(set, separated, replace = TRUE,
                                             bootstrap_samples = 2, seed = 2))
  clipped <- suppressWarnings(match_each_model(set, separated, replace = TRUE,
                                               bootstrap_samples = 2, seed = 2,
                                               clip_bootstrap = TRUE))
  treated <- separated$treat == 1
  expect_gt(length(unique(plain$bootstrap$rank[treated, 1])),
            length(unique(clipped$bootstrap$rank[treated, 1])))
  expect_true(clipped$settings$clip_bootstrap)
  expect_false(plain$settings$clip_bootstrap)
  # The default equals clip_bootstrap = FALSE given explicitly.
  explicit <- suppressWarnings(match_each_model(set, separated, replace = TRUE,
                                                bootstrap_samples = 2, seed = 2,
                                                clip_bootstrap = FALSE))
  expect_identical(explicit$bootstrap$rank, plain$bootstrap$rank)
})

test_that("glm warnings under separation are collected into one warning", {
  set.seed(5)
  n <- 300
  separated <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  separated$treat <- as.integer(separated$x1 > 0.5)
  run <- collect_warnings(fit_propensity_models(treat ~ x1 + x2, separated,
                                                algorithms = "logistic_regression"))
  expect_length(run$warnings, 1)
  expect_match(run$warnings, "glm\\(\\) warned in [0-9]+ of 5 cross-fitting folds")
  log <- run$value$logistic_log
  expect_true(any(log$glm_fitted_0_or_1 | log$glm_not_converged))
})

test_that("rank disagreement with equal weight per algorithm is the mean of d_between", {
  result <- quiet_gate(propensity_sensitivity(review_formula, data = lalonde,
                                              algorithms = review_fast, seeds = 1:2,
                                              n_trees = 50, bootstrap_samples = 0, seed = 3))
  table <- result$summary$table
  treated <- result$patients$treated == 1
  equal_weight <- table$with_gate[table$quantity == "mean_rank_disagreement_equal_weight_treated"]
  expect_equal(equal_weight, mean(result$patients$rank_disagreement_between[treated]))
  # By hand: the standard deviation over algorithms of each algorithm's mean rank.
  ranks <- result$components$model_set$rank[, result$summary$models_used, drop = FALSE]
  algorithms <- result$models$algorithm[match(colnames(ranks), result$models$model_id)]
  if (length(unique(algorithms)) >= 2) {
    means <- vapply(unique(algorithms), function(a) {
      rowMeans(ranks[, algorithms == a, drop = FALSE])
    }, numeric(nrow(ranks)))
    expect_equal(mean(apply(means, 1, stats::sd)[treated]), equal_weight)
  }
})

test_that("seed_sensitivity() reruns the diagnostic once per seed", {
  by_seed <- seed_sensitivity(review_formula, seeds = c(2, 5), data = lalonde,
                              algorithms = review_fast, model_seeds = 1, n_trees = 50,
                              bootstrap_samples = 2, keep_results = TRUE)
  expect_s3_class(by_seed, "seed_sensitivity")
  expect_identical(by_seed$table$seed, c(2, 5))
  expect_named(by_seed$table, c("seed", "models_kept", "models_fitted", "algorithms_kept",
                                "kept_by_algorithm", "fallback", "rank_disagreement",
                                "rank_disagreement_equal_weight", "pair_instability",
                                "bootstrap_rank_instability"))
  single <- quiet_gate(propensity_sensitivity(review_formula, data = lalonde,
                                              algorithms = review_fast, seeds = 1, n_trees = 50,
                                              bootstrap_samples = 2, seed = 5))
  table <- single$summary$table
  expect_equal(by_seed$table$rank_disagreement[2],
               table$with_gate[table$quantity == "mean_rank_disagreement_treated"])
  expect_equal(by_seed$range$minimum[1], min(by_seed$table$rank_disagreement))
  expect_output(print(by_seed), "Report every row")
  expect_error(seed_sensitivity(review_formula, data = lalonde, seed = 3), "do not pass `seed`")
})

test_that("bootstrap_rank_instability() with a formula equals the stored bootstrap ranks", {
  model_set <- fit_propensity_models(review_formula, lalonde, algorithms = "logistic_regression")
  stored <- match_each_model(model_set, lalonde, replace = TRUE, bootstrap_samples = 4, seed = 9)
  direct <- bootstrap_rank_instability(review_formula, data = lalonde, bootstrap_samples = 4,
                                       seed = 9)
  expect_equal(direct$summary$mean_treated,
               bootstrap_rank_instability(stored)$summary$mean_treated)
})

test_that("a factor made inside the formula is checked for unseen levels", {
  frame <- data.frame(x1 = seq(-1, 1, length.out = 40), code = c(9, rep(1:2, length.out = 39)))
  frame$treat <- rep(0:1, 20)
  prepared <- prepare_model_data(treat ~ x1 + factor(code), frame)
  train <- seq_len(40) != 1
  expect_true(has_unseen_levels(prepared, train, !train))
  expect_false(has_unseen_levels(prepare_model_data(treat ~ x1 + code, frame), train, !train))
})

test_that("a stopped parallel worker is reported with the task name", {
  skip_on_os("windows")
  # A task that returns NULL looks like a worker stopped by the system.
  expect_error(run_tasks(list(1, 2), function(task) if (task == 2) NULL else task, cores = 2,
                         task_names = c("first", "second"), what = "model fit"),
               "The model fit \"second\" returned no result")
})

test_that("the matching seeds are a family of their own", {
  model_seeds <- vapply(100L + 1:7, function(i) derive_seed(1, i), integer(1))
  matching_seeds <- vapply(matching_seed_offset + 1:7, function(i) derive_seed(1, i),
                           integer(1))
  bootstrap_draws <- vapply(1000L + 1:2000, function(i) derive_seed(1, i), integer(1))
  bootstrap_matching <- vapply(bootstrap_matching_seed_offset + 1:2000,
                               function(i) derive_seed(1, i), integer(1))
  expect_length(intersect(model_seeds, matching_seeds), 0)
  expect_length(intersect(bootstrap_draws, bootstrap_matching), 0)
})
