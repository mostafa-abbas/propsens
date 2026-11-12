data("lalonde", package = "MatchIt")
lalonde_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75
lalonde_covariates <- lalonde[, c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
treated_lalonde <- lalonde$treat == 1

# The analyst's own match, as an applied user would write it.
analyst_match <- function(replace = TRUE, ...) {
  # nolint next: object_usage_linter. `lalonde` is loaded by data() at the top of this file.
  MatchIt::matchit(lalonde_formula, data = lalonde, method = "nearest", distance = "glm",
                   link = "linear.logit", replace = replace, ...)
}

# Small settings keep the test suite fast: 5 models (logistic + 2 forests + 2 boosting). The
# wrapper's logic does not depend on the model set; the default set is tested in test-bart_grf.R.
fast_set <- c("logistic_regression", "random_forest", "gradient_boosting")
small_run <- function(x, ...) {
  sensitivity_quietly(x, algorithms = fast_set, seeds = 1:2, n_trees = 50,
                      bootstrap_samples = 5, ...)
}
with_replacement <- small_run(analyst_match(TRUE))
# The absolute rule, an option; on lalonde no model reaches 0.10.
with_replacement_absolute <- small_run(analyst_match(TRUE), balance_rule = "absolute")

test_that("1. a matchit object: formula and options are read, results equal the pieces", {
  settings <- with_replacement$settings
  expect_equal(settings$formula, lalonde_formula)
  expect_identical(settings$source, "matchit")
  expect_true(settings$replace)
  expect_identical(settings$m_order, "largest")
  expect_null(settings$caliper_sd)
  expect_identical(settings$covariate_names, names(lalonde_covariates))

  smallest <- read_matchit_design(analyst_match(FALSE, m.order = "smallest"), lalonde,
                                  environment())
  expect_identical(smallest$m_order, "smallest")
  expect_false(smallest$replace)

  model_set <- fit_propensity_models(lalonde_formula, lalonde, seeds = 1:2, n_trees = 50,
                                     seed = 1, algorithms = fast_set)
  match_set <- match_each_model(model_set, lalonde, replace = TRUE, m_order = "largest",
                                bootstrap_samples = 5, seed = 1)
  balance <- balance_spread(match_set, model_set$treatment, lalonde_covariates,
                            use_gate = FALSE)
  # Default relative rule: within 0.10 of the best model's max |SMD|, and <= 0.25.
  smd <- balance$models$smd_max_main
  passed <- model_set$models$passed_logloss_gate & smd <= min(smd) + 0.10 & smd <= 0.25
  ids <- if (sum(passed) < 2) model_set$models$model_id else model_set$models$model_id[passed]
  by_hand_rank <- rank_disagreement(model_set, models = ids)
  by_hand_control <- control_spread(match_set, lalonde_covariates, models = ids)
  by_hand_recurrence <- pair_recurrence(match_set, models = ids)

  patients <- with_replacement$patients
  expect_identical(with_replacement$summary$models_used, ids)
  expect_equal(patients$rank_disagreement, by_hand_rank$patients$rank_disagreement)
  expect_equal(patients$rank_disagreement_within,
               by_hand_rank$patients$rank_disagreement_within)
  expect_equal(patients$control_spread_scaled, by_hand_control$patients$control_spread_scaled)
  expect_equal(patients$pair_recurrence_models,
               by_hand_recurrence$patients$pair_recurrence_models)
  expect_equal(patients$pair_recurrence_bootstrap,
               by_hand_recurrence$patients$pair_recurrence_bootstrap)
  models <- with_replacement$models[!with_replacement$models$is_reference, ]
  expect_equal(models$smd_max_main, balance$models$smd_max_main)
  expect_equal(models$logloss, model_set$models$logloss)
  table <- with_replacement$summary$table
  expect_equal(table$with_gate[table$quantity == "mean_rank_disagreement_treated"],
               by_hand_rank$summary$mean_d_treated)
  expect_equal(table$with_gate[table$quantity == "mean_control_spread_scaled_treated"],
               by_hand_control$summary$mean_scaled_treated)
})

test_that("2. without replacement; with a caliper the treated-kept share is the object's", {
  without <- small_run(analyst_match(FALSE))
  expect_false(without$settings$replace)
  expect_false(without$components$match_set$settings$replace)
  expect_true(all(without$models$treated_kept == 1))

  analyst <- analyst_match(FALSE, caliper = 0.1)
  object_share <- sum(analyst$weights[treated_lalonde] > 0) / sum(treated_lalonde)
  expect_lt(object_share, 1)
  # cross_fit = FALSE: the logistic model is the analyst's in-sample model, so it must give
  # exactly the analyst's match (the caliper is converted to standard deviations).
  in_sample <- sensitivity_quietly(analyst, algorithms = fast_set, cross_fit = FALSE, seeds = 1,
                                   n_trees = 50,
                                   bootstrap_samples = 0)
  expect_equal(in_sample$settings$caliper_sd, 0.1)
  models <- in_sample$models
  expect_equal(models$treated_kept[models$model_id == "reference"], object_share)
  expect_equal(models$treated_kept[models$model_id == "logistic_regression_1"], object_share)
  ours <- in_sample$components$match_set$matches$logistic_regression_1
  expect_identical(ours$control_id,
                   unname(analyst$match.matrix[!is.na(analyst$match.matrix[, 1]), 1]))
  table <- in_sample$summary$table
  expect_equal(table$with_gate[table$quantity == "treated_kept_share_min"],
               min(models$treated_kept[models$in_use]))
})

test_that("3. clear errors: other methods, ratio 2, missing data, data not found", {
  full_matching <- analyst_match(TRUE)
  full_matching$info$method <- "full"
  expect_error(propensity_sensitivity(full_matching),
               "only 1:1 nearest-neighbour matching is supported in this version")
  ratio_two <- analyst_match(TRUE, ratio = 2)
  expect_error(propensity_sensitivity(ratio_two),
               "only 1:1 nearest-neighbour matching is supported in this version")
  with_missing <- lalonde
  with_missing$age[3] <- NA
  expect_error(propensity_sensitivity(lalonde_formula, with_missing), "Missing values")
  make_object_with_local_data <- function() {
    local_data <- lalonde
    MatchIt::matchit(lalonde_formula, data = local_data, distance = "glm",
                     link = "linear.logit", replace = TRUE)
  }
  expect_error(propensity_sensitivity(make_object_with_local_data()), "pass the data")
  expect_error(propensity_sensitivity(analyst_match(TRUE), replace = FALSE),
               "read from the matchit object")
  expect_error(propensity_sensitivity(lalonde_formula), "pass the data")
  expect_error(propensity_sensitivity("treat"), "matchit object .* or a formula")
})

test_that("4. the reference row: log-loss NA, never in use, flagged when it fails balance", {
  with_replacement <- with_replacement_absolute
  models <- with_replacement$models
  reference <- models[models$model_id == "reference", ]
  expect_equal(nrow(reference), 1)
  expect_true(reference$is_reference)
  expect_false(reference$in_use)
  expect_true(is.na(reference$logloss))
  expect_true(is.na(reference$passed_logloss_gate))
  expect_false("reference" %in% with_replacement$summary$models_used)
  analyst <- analyst_match(TRUE)
  balance <- balance_statistics(lalonde_covariates, lalonde$treat,
                                list(reference = as.numeric(analyst$weights)))
  expect_equal(reference$smd_max_main, balance$smd_max_main)
  # On lalonde the analyst's match fails the 0.10 balance gate: a real failing case.
  expect_gt(reference$smd_max_main, 0.10)
  expect_false(reference$passed_quality_gate)
  expect_match(reference$gate_note, "judged on balance only: above the balance limit",
               fixed = TRUE)
  expect_true(with_replacement$summary$reference_model_not_kept)
  expect_match(printed_text(with_replacement),
               "Your own model (glm, linear.logit): above the balance limit", fixed = TRUE)
  passing <- reference_model_row(
    list(label = "glm, linear.logit", rank = rep(0.5, nrow(lalonde)),
         weights = as.numeric(analyst$weights)),
    lalonde_covariates, lalonde$treat, balance_limit = 1, columns = names(models)
  )
  expect_true(passing$passed_quality_gate)
  expect_match(passing$gate_note, "judged on balance only: within the balance limit",
               fixed = TRUE)

  reference_rank <- rank(as.numeric(analyst$distance)) / nrow(lalonde)
  ranks <- with_replacement$components$model_set$rank
  ids <- with_replacement$summary$models_used
  table <- with_replacement$summary$table
  expect_equal(table$with_gate[table$quantity ==
                                 "reference_mean_absolute_rank_difference_treated"],
               mean(abs(ranks[treated_lalonde, ids] - reference_rank[treated_lalonde])))
})

test_that("5. the final gate uses both parts; the fewer-than-two rule works on the final set", {
  log_loss_table <- data.frame(model_id = c("a", "b", "c", "d"),
                               passed_logloss_gate = c(TRUE, TRUE, TRUE, FALSE))
  balance_table <- data.frame(model_id = c("a", "b", "c", "d"),
                              smd_max_main = c(0.05, 0.20, 0.08, 0.01),
                              smd_max_extended = 0, energy_distance = 0, ks_max = 0)
  absolute <- function(threshold) list(rule = "absolute", threshold = threshold)
  gate <- apply_final_quality_gate(log_loss_table, balance_table, absolute(0.10))
  expect_identical(gate$final_ids, c("a", "c"))
  expect_false(gate$fewer_than_two)
  expect_match(gate$models$gate_note[2], "max \\|SMD\\| main terms above balance_threshold")
  expect_match(gate$models$gate_note[4], "log-loss above")
  single <- apply_final_quality_gate(log_loss_table, balance_table, absolute(0.06))
  expect_true(single$fewer_than_two)
  expect_identical(single$final_ids, c("a", "b", "c", "d"))

  # On lalonde: every model passes log-loss (tolerance 1); the worst-balanced one is excluded.
  models <- with_replacement$models[!with_replacement$models$is_reference, ]
  threshold <- sort(models$smd_max_main)[nrow(models) - 1]
  worst <- models$model_id[which.max(models$smd_max_main)]
  gated <- small_run(analyst_match(TRUE), quality_tolerance = 1, balance_rule = "absolute",
                     balance_threshold = threshold)
  gated_models <- gated$models
  expect_true(gated_models$passed_logloss_gate[gated_models$model_id == worst])
  expect_false(gated_models$passed_quality_gate[gated_models$model_id == worst])
  expect_false(worst %in% gated$summary$models_used)
  expect_length(gated$summary$models_used, nrow(models) - 1)
  expect_false(gated$quality_gate_left_fewer_than_two_models)
  expect_equal(gated$patients$rank_disagreement,
               rank_disagreement(gated$components$model_set,
                                 models = gated$summary$models_used)$patients$rank_disagreement)
  # With the absolute 0.10 rule no lalonde model passes balance: all are used, flag set.
  expect_true(with_replacement_absolute$quality_gate_left_fewer_than_two_models)
  expect_identical(with_replacement_absolute$summary$models_used, models$model_id)
})

test_that("5b. the relative balance rule: margin above the best model, and a cap", {
  log_loss_table <- data.frame(model_id = c("a", "b", "c", "d", "e"),
                               passed_logloss_gate = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  balance_table <- data.frame(model_id = c("a", "b", "c", "d", "e"),
                              smd_max_main = c(0.12, 0.21, 0.23, 0.30, 0.01),
                              smd_max_extended = 0, energy_distance = 0, ks_max = 0)
  relative <- function(margin = 0.10, cap = 0.25) {
    list(rule = "relative", margin = margin, cap = cap)
  }
  # The best model is e (0.01; it fails log-loss but still sets the best balance): limit
  # 0.11, so none of a-d passes balance and the fewer-than-two rule applies.
  gate <- apply_final_quality_gate(log_loss_table, balance_table, relative())
  expect_equal(gate$balance_limit, 0.11)
  expect_true(gate$fewer_than_two)
  # Without e the best is a (0.12): limit 0.22 -> a and b pass; c (0.23) is more than 0.10 above.
  without_e <- apply_final_quality_gate(log_loss_table[1:4, ], balance_table[1:4, ], relative())
  expect_equal(without_e$balance_limit, 0.22)
  expect_identical(without_e$final_ids, c("a", "b"))
  expect_match(without_e$models$gate_note[3], "more than balance_margin above the best")
  # The cap binds: margin 0.5 would allow 0.62, the cap 0.25 excludes d (0.30).
  capped <- apply_final_quality_gate(log_loss_table[1:4, ], balance_table[1:4, ],
                                     relative(margin = 0.5))
  expect_equal(capped$balance_limit, 0.25)
  expect_identical(capped$final_ids, c("a", "b", "c"))
  # Equal to the simulation's rule (simulation/run_one_dataset.R, final_gate_sets).
  smd <- balance_table$smd_max_main[1:4]
  expect_identical(without_e$models$passed_balance_gate, smd <= min(smd) + 0.10 & smd <= 0.25)

  # On lalonde: default rule recorded in settings, the limit equals the best model + 0.10 or the
  # cap, and print() and summary() state the rule in effect; the absolute rule says so too.
  settings <- with_replacement$settings
  expect_identical(settings$balance_rule, "relative")
  fitted <- with_replacement$models[!with_replacement$models$is_reference, ]
  expect_equal(settings$balance_limit, min(min(fitted$smd_max_main) + 0.10, 0.25))
  expect_identical(fitted$passed_balance_gate, fitted$smd_max_main <= settings$balance_limit)
  expect_output(print(with_replacement), "relative balance rule")
  expect_output(print(summary(with_replacement)), "relative balance rule")
  expect_match(printed_text(with_replacement_absolute),
               "absolute balance rule: max |SMD| main terms at most 0.10", fixed = TRUE)
  expect_output(print(summary(with_replacement_absolute)), "absolute balance rule")
  expect_error(small_run(analyst_match(TRUE), balance_margin = -1), "balance_margin")

  # "relative_only" drops the cap: margin 0.5 would allow 0.62, so d (0.30) now passes.
  uncapped <- apply_final_quality_gate(log_loss_table[1:4, ], balance_table[1:4, ],
                                       list(rule = "relative_only", margin = 0.5, cap = 0.25))
  expect_equal(uncapped$balance_limit, 0.62)
  expect_identical(uncapped$final_ids, c("a", "b", "c", "d"))
  # With the default margin the cap does not bind here, so both rules agree.
  relative_only_default <- apply_final_quality_gate(log_loss_table[1:4, ], balance_table[1:4, ],
                                                    list(rule = "relative_only", margin = 0.10))
  expect_identical(relative_only_default$final_ids, without_e$final_ids)
  expect_match(uncapped$models$gate_note[1], "passed")
  expect_error(small_run(analyst_match(TRUE), balance_rule = "loose"), "should be one of")
  relative_only <- small_run(analyst_match(TRUE), balance_rule = "relative_only")
  fitted_only <- relative_only$models[!relative_only$models$is_reference, ]
  expect_equal(relative_only$settings$balance_limit, min(fitted_only$smd_max_main) + 0.10)
  expect_output(print(relative_only), "relative-only balance rule \\(no cap\\)")
})

test_that("6. cross_fit = FALSE works and gives in-sample logistic scores", {
  analyst <- analyst_match(TRUE)
  in_sample <- sensitivity_quietly(analyst, algorithms = fast_set, cross_fit = FALSE, seeds = 1,
                                   n_trees = 50,
                                   bootstrap_samples = 0)
  expect_false(in_sample$settings$cross_fit)
  expect_equal(unname(in_sample$patients$logistic_score), unname(analyst$distance),
               tolerance = 1e-6)
  expect_true(all(is.na(in_sample$patients$pair_recurrence_bootstrap)))
})

test_that("7. the formula interface equals the matchit-object interface", {
  by_formula <- small_run(lalonde_formula, data = lalonde, replace = TRUE)
  expect_identical(by_formula$settings$source, "formula")
  expect_equal(by_formula$patients, with_replacement$patients)
  object_models <- with_replacement$models[!with_replacement$models$is_reference, ]
  expect_equal(by_formula$models, object_models)
  expect_false(any(by_formula$models$is_reference))
  rows <- by_formula$summary$table$quantity != "reference_mean_absolute_rank_difference_treated"
  expect_equal(by_formula$summary$table[rows, ], with_replacement$summary$table[rows, ])
})

test_that("8. print, summary and as.data.frame", {
  # The printed snapshot uses an analyst's match on the glm logit rounded to 8 digits. lalonde
  # has controls with identical covariates, hence exactly tied scores; last-digit differences
  # in the glm fit between BLAS libraries (reference BLAS, OpenBLAS, macOS) would otherwise
  # decide which tied control the analyst's match uses. Rounded, the ties are exact and
  # MatchIt breaks them the same way on every platform.
  rounded_logit <- round(stats::predict(stats::glm(lalonde_formula, data = lalonde,
                                                   family = stats::binomial())), 8)
  rounded_match <- MatchIt::matchit(lalonde_formula, data = lalonde, method = "nearest",
                                    distance = rounded_logit, replace = TRUE)
  rounded_run <- small_run(rounded_match)
  expect_snapshot(print(rounded_run))
  expect_snapshot(print(summary(rounded_run)))
  as_frame <- as.data.frame(with_replacement)
  expect_equal(as_frame, with_replacement$patients, ignore_attr = "note")
  expect_match(attr(as_frame, "note"), "not flags for individual patients")
  expect_named(with_replacement$patients,
               c("id", "treated", "rank_disagreement", "rank_disagreement_between",
                 "rank_disagreement_within", "logit_disagreement", "control_spread",
                 "control_spread_scaled", "pair_recurrence_models",
                 "pair_recurrence_bootstrap", "outside_support", "logistic_score"))
  controls <- with_replacement$patients[with_replacement$patients$treated == 0, ]
  expect_true(all(is.na(controls$control_spread_scaled)))
  expect_true(all(is.na(controls$outside_support)))
})

test_that("the four diagnostics accept a set of models", {
  model_set <- with_replacement$components$model_set
  match_set <- with_replacement$components$match_set
  ids <- c("logistic_regression_1", "random_forest_1", "random_forest_2")
  chosen <- rank_disagreement(model_set, models = ids)
  expect_identical(chosen$models_used, ids)
  expect_true(chosen$settings$models_supplied)
  expect_identical(control_spread(match_set, lalonde_covariates, models = ids)$models_used, ids)
  expect_identical(pair_recurrence(match_set, models = ids)$models_used, ids)
  balance <- balance_spread(match_set, model_set$treatment, lalonde_covariates, models = ids)
  expect_identical(balance$models_used, ids)
  expect_identical(balance$models$in_use, balance$models$model_id %in% ids)
  expect_error(rank_disagreement(model_set, models = "no_such_model"), "model ids")
  expect_true(any(grepl("models chosen by the caller", capture.output(print(chosen)))))
})

test_that("9. timing: lalonde, package defaults, one core, under 60 seconds", {
  skip_unless_timing_tests()
  analyst <- analyst_match(TRUE)
  timing <- system.time(result <- sensitivity_quietly(analyst, cores = 1))
  expect_elapsed_under(timing, 60, "propensity_sensitivity, lalonde, defaults")
  expect_s3_class(result, "propensity_sensitivity")
  expect_equal(nrow(result$models), 8)
})

test_that("print and summary lead with gate composition, rank disagreement, balance spread", {
  printed <- utils::capture.output(print(with_replacement))
  summarised <- utils::capture.output(print(summary(with_replacement)))
  expect_identical(summarised[1], "Propensity model sensitivity: summary")
  for (lines in list(printed, summarised)) {
    gate <- grep("^1\\. Quality gate:", lines)
    rank <- grep("^2\\. Rank disagreement", lines)
    spread <- grep("^3\\. Balance spread", lines)
    expect_length(gate, 1)
    expect_true(gate < rank && rank < spread)
    # The composition names the number of algorithms kept and lists the excluded models.
    used_ids <- with_replacement$summary$models_used
    in_use <- with_replacement$models[with_replacement$models$model_id %in% used_ids, ]
    expect_match(lines[gate], paste0(length(unique(in_use$algorithm)), " algorithm"))
    excluded <- setdiff(with_replacement$models$model_id[!with_replacement$models$is_reference],
                        with_replacement$summary$models_used)
    for (id in excluded) expect_true(any(grepl(paste0("^     ", id, ": "), lines)))
    # Every printed line fits the console width.
    expect_true(all(nchar(lines[seq_len(spread)]) <= getOption("width")))
  }
  table <- with_replacement$summary$table
  joined <- printed_text(with_replacement)
  expect_match(joined, paste0("3. Balance spread (descriptive; range of max |SMD| main terms ",
                              "over the kept models): ",
                              format_number(table$with_gate[table$quantity ==
                                                              "smd_max_main_range"], 3)),
               fixed = TRUE)
  # The share of treated kept by each kept model is shown.
  for (id in with_replacement$summary$models_used) {
    expect_match(joined, paste0(id, "\\s+100 %"))
  }
})

test_that("the fitted model set is returned, equal to fit_propensity_models()", {
  # Callers (for example, other packages) reuse the scores without refitting.
  sensitivity <- sensitivity_quietly(treat ~ age + educ + race + married + nodegree + re74 +
                                       re75, data = lalonde, seeds = 1, n_trees = 50,
                                     algorithms = fast_set, bootstrap_samples = 2,
                                     replace = FALSE, seed = 3)
  direct <- fit_propensity_models(treat ~ age + educ + race + married + nodegree + re74 + re75,
                                  lalonde, seeds = 1, n_trees = 50, seed = 3,
                                  algorithms = fast_set)
  returned <- sensitivity$components$model_set
  expect_s3_class(returned, "propensity_model_set")
  expect_identical(returned$probability, direct$probability)
  expect_identical(returned$models$model_id, direct$models$model_id)
  expect_true(all(sensitivity$summary$models_used %in% colnames(returned$probability)))
})

test_that("the gate fallback is shown first, in plain words", {
  expect_true(with_replacement_absolute$quality_gate_left_fewer_than_two_models)
  printed <- utils::capture.output(print(with_replacement_absolute))
  summarised <- utils::capture.output(print(summary(with_replacement_absolute)))
  expect_match(printed[2], "^WARNING: fewer than two models passed the quality gate, so there")
  expect_match(summarised[2], "^WARNING: fewer than two models passed the quality gate")
  fallback_words <- paste("there is no set of similarly good models to compare. Every model is",
                          "used below; read the values as a description of all models, not of",
                          "a kept set.")
  expect_match(printed_text(with_replacement_absolute), fallback_words, fixed = TRUE)
  expect_lt(grep("^WARNING", printed), grep("^1\\. Quality gate:", printed))
  expect_false(with_replacement$quality_gate_left_fewer_than_two_models)
  expect_false(any(grepl("WARNING", utils::capture.output(print(with_replacement)))))
  expect_false(any(grepl("WARNING", utils::capture.output(print(summary(with_replacement))))))
  # With the fallback the summary prints one column (every model is used).
  expect_true(any(grepl("every model is used, so one column", summarised, fixed = TRUE)))
})

test_that("bootstrap rank instability is in the output, next to rank disagreement", {
  table <- with_replacement$summary$table
  row <- table$quantity == "bootstrap_rank_instability_treated"
  value <- with_replacement$bootstrap_rank_instability$summary$mean_treated
  expect_equal(table$with_gate[row], value)
  expect_equal(table$without_gate[row], value)
  ranks <- with_replacement$components$match_set$bootstrap$rank
  expect_equal(value, mean(apply(ranks[treated_lalonde, ], 1, stats::sd)))
  for (lines in list(utils::capture.output(print(with_replacement)),
                     utils::capture.output(print(summary(with_replacement))))) {
    rank <- grep("^2\\. Rank disagreement", lines)
    spread <- grep("^3\\. Balance spread", lines)
    bootstrap <- grep("Bootstrap rank instability of the logistic model \\(descriptive\\)", lines)
    expect_true(rank < bootstrap && bootstrap < spread)
    expect_match(lines[bootstrap], format_number(value, 3), fixed = TRUE)
  }
})

test_that("pair instability = 1 - pair recurrence, shown with its maximum", {
  table <- with_replacement$summary$table
  instability <- table[table$quantity == "pair_instability_models_mean_treated", ]
  recurrence <- table[table$quantity == "pair_recurrence_models_mean_treated", ]
  expect_equal(nrow(instability), 1)
  expect_equal(nrow(recurrence), 1)
  expect_equal(instability$with_gate, 1 - recurrence$with_gate)
  expect_equal(instability$without_gate, 1 - recurrence$without_gate)
  patients <- with_replacement$patients
  expect_equal(instability$with_gate,
               1 - mean(patients$pair_recurrence_models[patients$treated == 1]))
  n_kept <- length(with_replacement$summary$models_used)
  expect_equal(with_replacement$summary$pair_instability_maximum, 1 - 1 / n_kept)
  for (lines in list(utils::capture.output(print(with_replacement)),
                     utils::capture.output(print(summary(with_replacement))))) {
    joined <- gsub("\\s+", " ", paste(lines, collapse = " "))
    expect_match(joined, paste0("Pair instability (exploratory; study mean over treated): ",
                                format_number(instability$with_gate, 3), " of a possible ",
                                format_number(1 - 1 / n_kept, 3), " with ", n_kept, " models"),
                 fixed = TRUE)
    expect_lt(grep("Equal weight per algorithm", lines), grep("^   Pair instability", lines))
  }
})

test_that("without bootstrap refits the instability is NA and printing still works", {
  none <- sensitivity_quietly(lalonde_formula, lalonde, algorithms = fast_set, seeds = 1,
                              n_trees = 50, bootstrap_samples = 0, replace = FALSE)
  expect_null(none$bootstrap_rank_instability)
  table <- none$summary$table
  expect_true(is.na(table$with_gate[table$quantity == "bootstrap_rank_instability_treated"]))
  expect_true(any(grepl("instability of the logistic model.*NA",
                        utils::capture.output(print(none)))))
})

test_that("with a formula the default is 1:1 matching without replacement", {
  expect_false(eval(formals(propensity_sensitivity)$replace))
  by_default <- sensitivity_quietly(lalonde_formula, lalonde, algorithms = fast_set, seeds = 1,
                                    n_trees = 50, bootstrap_samples = 0, seed = 5)
  expect_false(by_default$settings$replace)
  explicit <- sensitivity_quietly(lalonde_formula, lalonde, algorithms = fast_set, seeds = 1,
                                  n_trees = 50, bootstrap_samples = 0, replace = FALSE, seed = 5)
  expect_equal(by_default$patients, explicit$patients)
  # A matchit object keeps its own option.
  expect_true(with_replacement$settings$replace)
})
