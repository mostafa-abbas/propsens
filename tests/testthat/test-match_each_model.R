data("lalonde", package = "MatchIt")
lalonde_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75
# "As MatchIt does it": one in-sample logistic model, so that MatchIt's own glm match is the
# reference (test 1). The second set is a small cross-fitted set with all three algorithms.
logistic_set <- fit_propensity_models(lalonde_formula, lalonde, algorithms = "logistic_regression",
                                      cross_fit = FALSE)
model_set <- fit_propensity_models(lalonde_formula, lalonde, seeds = 1, n_trees = 50,
                                   algorithms = faster_set)
n_lalonde <- nrow(lalonde)
treated_lalonde <- lalonde$treat == 1

reference_matchit <- function(replace, caliper = NULL) {
  # nolint next: object_usage_linter. `lalonde` is loaded by data() at the top of this file.
  MatchIt::matchit(lalonde_formula, data = lalonde, method = "nearest", distance = "glm",
                   link = "linear.logit", replace = replace, m.order = "largest",
                   caliper = caliper)
}

# Pairs that differ from MatchIt's are allowed only as near-tie swaps: the treated patient is
# equally far, within `tolerance` on MatchIt's logit, from our partner and from MatchIt's.
# Such ties are broken by the last bits of floating-point arithmetic, which differ between
# compilers (for example NSW71, halfway between PSID392 and PSID400 in lalonde).
near_tie_swaps <- function(control_id, reference, tolerance = 1e-8) {
  # nolint next: object_usage_linter. `lalonde` is loaded by data() at the top of this file.
  logit <- stats::setNames(as.numeric(reference$distance), rownames(lalonde))
  treated_id <- rownames(reference$match.matrix)
  expected <- unname(reference$match.matrix[, 1])
  differ <- which(is.na(control_id) != is.na(expected) |
                    (!is.na(control_id) & !is.na(expected) & control_id != expected))
  gap <- abs(abs(logit[treated_id[differ]] - logit[control_id[differ]]) -
               abs(logit[treated_id[differ]] - logit[expected[differ]]))
  list(rows = differ, all_near_ties = !anyNA(gap) && all(gap <= tolerance),
       controls = unique(c(control_id[differ], expected[differ])))
}

test_that("1. reproduction: the same pairs and weights as MatchIt's own glm match", {
  for (replace in c(FALSE, TRUE)) {
    ours <- match_each_model(logistic_set, lalonde, replace = replace)
    reference <- reference_matchit(replace)
    pairs <- ours$matches$logistic_regression_1
    expect_identical(pairs$treated_id, rownames(reference$match.matrix))
    # Identical pairs, except near-tie swaps, and at most two of them.
    swaps <- near_tie_swaps(pairs$control_id, reference)
    expect_true(swaps$all_near_ties)
    expect_lte(length(swaps$rows), 2)
    # Weights identical, except for the controls of a near-tie swap.
    unaffected <- !rownames(lalonde) %in% swaps$controls
    if (!replace) {
      expect_equal(unname(ours$weights[unaffected, 1]), unname(reference$weights[unaffected]))
    } else {
      # MatchIt normalizes control weights to a mean of 1 by default; ours are the raw counts
      # (normalize = FALSE), so the control weights are proportional.
      kept <- !treated_lalonde & unaffected
      ours_control <- ours$weights[kept, 1]
      reference_control <- reference$weights[kept]
      expect_equal(unname(ours_control * sum(reference_control) / sum(ours_control)),
                   unname(reference_control))
      expect_equal(unname(ours$weights[treated_lalonde, 1]),
                   unname(reference$weights[treated_lalonde]))
    }
  }
})

test_that("1b. the near-tie rule accepts the NSW71 swap and rejects any other change", {
  reference <- reference_matchit(TRUE)
  control_id <- unname(reference$match.matrix[, 1])
  expect_identical(near_tie_swaps(control_id, reference)$rows, integer(0))
  row <- which(rownames(reference$match.matrix) == "NSW71")
  other <- setdiff(c("PSID392", "PSID400"), control_id[row])
  tie_swapped <- replace(control_id, row, other)
  expect_true(near_tie_swaps(tie_swapped, reference)$all_near_ties)
  expect_setequal(near_tie_swaps(tie_swapped, reference)$controls, c("PSID392", "PSID400"))
  far_swapped <- replace(control_id, row, "PSID1")
  expect_false(near_tie_swaps(far_swapped, reference)$all_near_ties)
})

test_that("2. without replacement: every treated matched once; a caliper drops treated", {
  ours <- match_each_model(model_set, lalonde, replace = FALSE)
  n_treated <- sum(treated_lalonde)
  for (id in ours$model_ids) {
    pairs <- ours$matches[[id]]
    expect_equal(nrow(pairs), n_treated)
    expect_equal(sort(pairs$treated_id), sort(rownames(lalonde)[treated_lalonde]))
    expect_false(any(duplicated(pairs$control_id)))
    expect_equal(unname(ours$treated_kept[id]), 1)
    expect_true(all(ours$weights[, id] %in% c(0, 1)))
  }
  with_caliper <- match_each_model(logistic_set, lalonde, replace = FALSE, caliper = 0.05)
  reference <- reference_matchit(FALSE, caliper = 0.05)
  counts <- summary(reference)$nn
  expect_lt(unname(with_caliper$treated_kept[1]), 1)
  expect_equal(unname(with_caliper$treated_kept[1]),
               counts["Matched", "Treated"] / counts["All", "Treated"])
  expect_equal(nrow(with_caliper$matches[[1]]), counts["Matched", "Treated"])
  expect_identical(with_caliper$matches[[1]]$control_id,
                   unname(reference$match.matrix[!is.na(reference$match.matrix[, 1]), 1]))
})

test_that("3. with replacement: controls repeat, control weights sum to the matched treated", {
  ours <- match_each_model(model_set, lalonde, replace = TRUE)
  for (id in ours$model_ids) {
    pairs <- ours$matches[[id]]
    expect_true(any(duplicated(pairs$control_id)))
    expect_equal(sum(ours$weights[!treated_lalonde, id]), nrow(pairs))
    expect_equal(sum(ours$weights[!treated_lalonde, id]), sum(treated_lalonde))
    expect_true(all(ours$weights[treated_lalonde, id] == 1))
    counts <- table(pairs$control_id)
    expect_equal(unname(ours$weights[names(counts), id]), as.numeric(counts))
    expect_true(all(pairs$weight == 1))
    expect_equal(pairs$pair_id, seq_len(nrow(pairs)))
  }
})

test_that("4. oracle: present only when oracle_logit is given, matched with the same options", {
  without <- match_each_model(model_set, lalonde, replace = TRUE)
  expect_null(without$oracle)
  expect_false(without$settings$has_oracle)
  oracle_logit <- logistic_set$logit[, "logistic_regression_1"]
  with_oracle <- match_each_model(model_set, lalonde, replace = TRUE, oracle_logit = oracle_logit)
  expect_true(with_oracle$settings$has_oracle)
  expect_named(with_oracle$oracle, c("matches", "weights", "treated_kept", "call"))
  expect_false("oracle" %in% colnames(with_oracle$weights))
  # The oracle logit here is the in-sample logistic logit, so the oracle matches must equal
  # the logistic matches of the "as MatchIt does it" set under the same options.
  logistic_matches <- match_each_model(logistic_set, lalonde, replace = TRUE)
  expect_identical(with_oracle$oracle$matches, logistic_matches$matches$logistic_regression_1)
  expect_equal(with_oracle$oracle$weights, logistic_matches$weights[, 1])
  options_of <- function(text) sub(".*, replace = ", "replace = ", text)
  expect_identical(options_of(with_oracle$oracle$call),
                   options_of(with_oracle$calls[["logistic_regression_1"]]))
  expect_match(with_oracle$oracle$call, "distance = oracle_logit")
  expect_output(print(with_oracle), "oracle")
})

test_that("5. bootstrap: shapes, ranks, presence, original ids, B = 0, reproducibility", {
  n_samples <- 5
  ours <- match_each_model(logistic_set, lalonde, replace = TRUE, bootstrap_samples = n_samples,
                           seed = 7)
  boot <- ours$bootstrap
  expect_named(boot, c("rank", "pairs", "present", "skipped", "design_rule_samples",
                       "prediction_fallback_samples"))
  expect_equal(dim(boot$rank), c(n_lalonde, n_samples))
  expect_equal(dim(boot$pairs), c(n_lalonde, n_samples))
  expect_equal(dim(boot$present), c(n_lalonde, n_samples))
  expect_equal(nrow(boot$skipped), 0)
  expect_true(all(boot$rank >= 1 / n_lalonde & boot$rank <= 1))
  expect_identical(rownames(boot$rank), rownames(lalonde))
  for (b in seq_len(n_samples)) {
    rows <- propsens:::with_task_seed(propsens:::derive_seed(7, 1000 + b),
                                      sample.int(n_lalonde, n_lalonde, replace = TRUE))
    expected_present <- logical(n_lalonde)
    expected_present[rows] <- TRUE
    expect_identical(unname(boot$present[, b]), expected_present)
    # The refit on this resample reproduces the stored ranks for the original patients.
    refit <- glm(lalonde_formula, data = lalonde[rows, ], family = binomial())
    logit <- predict(refit, newdata = lalonde, type = "link")
    expect_equal(unname(boot$rank[, b]), unname(rank(logit, ties.method = "average")) / n_lalonde)
    pairs <- boot$pairs[, b]
    expect_true(all(is.na(pairs[!treated_lalonde])))
    expect_true(all(is.na(pairs[!boot$present[, b]])))
    expect_true(all(pairs[treated_lalonde & boot$present[, b]] %in% rownames(lalonde)))
    expect_true(all(lalonde[pairs[!is.na(pairs)], "treat"] == 0))
  }
  expect_null(match_each_model(logistic_set, lalonde, replace = TRUE)$bootstrap)
  again <- match_each_model(logistic_set, lalonde, replace = TRUE, bootstrap_samples = n_samples,
                            seed = 7)
  expect_identical(again$bootstrap, boot)
  expect_identical(again$matches, ours$matches)
  other_seed <- match_each_model(logistic_set, lalonde, replace = TRUE,
                                 bootstrap_samples = n_samples, seed = 8)
  expect_false(identical(other_seed$bootstrap$present, boot$present))
})

test_that("5b. bootstrap: a resample with fewer than 10 treated is skipped and logged", {
  # 10 treated out of 100: about half of the resamples have fewer than 10 treated.
  small <- data.frame(treat = rep(c(1, 0), c(10, 90)), x1 = rnorm(100), x2 = rnorm(100))
  small_set <- fit_propensity_models(treat ~ x1 + x2, small, algorithms = "logistic_regression",
                                     folds = 2)
  ours <- match_each_model(small_set, small, replace = TRUE, bootstrap_samples = 12, seed = 3)
  skipped <- ours$bootstrap$skipped
  expect_gt(nrow(skipped), 0)
  expect_match(skipped$reason, "fewer than 10")
  expect_equal(dim(ours$bootstrap$rank), c(100, 12))
  for (b in skipped$sample) {
    expect_true(all(is.na(ours$bootstrap$rank[, b])))
    expect_true(all(is.na(ours$bootstrap$present[, b])))
  }
  kept <- setdiff(seq_len(12), skipped$sample)
  expect_false(anyNA(ours$bootstrap$rank[, kept]))
  expect_output(print(ours), "skipped")
})

test_that("5c. bootstrap: a level absent from the resample uses the design-matrix rule", {
  rare <- lalonde
  # One patient with a level of its own: about 37 % of the resamples miss it.
  rare$region <- factor(c("a", rep(c("b", "c"), length.out = nrow(rare) - 1)))
  rare_formula <- update(lalonde_formula, . ~ . + region)
  rare_set <- fit_propensity_models(rare_formula, rare, algorithms = "logistic_regression",
                                    cross_fit = FALSE)
  ours <- match_each_model(rare_set, rare, replace = TRUE, bootstrap_samples = 8, seed = 5)
  boot <- ours$bootstrap
  expect_equal(nrow(boot$skipped), 0)
  expect_gt(length(boot$design_rule_samples), 0)
  expect_false(anyNA(boot$rank))
  # Per sample: the design rule is used exactly when the level is absent, and then equals a
  # direct glm.fit on the full-data design with the coefficient of the absent level set to 0;
  # a sample that has the level is refitted with glm() on the formula.
  prepared <- prepare_model_data(rare_formula, rare)
  n <- nrow(rare)
  for (b in seq_len(8)) {
    rows <- with_task_seed(derive_seed(5, 1000L + b), sample.int(n, n, replace = TRUE))
    resample <- rare[rows, , drop = FALSE]
    rownames(resample) <- NULL
    new <- bootstrap_sample_logits(rare_formula, rare, resample, rows, prepared)
    absent <- !1L %in% rows
    expect_identical(new$design_rule, absent)
    expect_identical(b %in% boot$design_rule_samples, absent)
    if (absent) {
      design <- cbind(1, prepared$design)
      direct <- stats::glm.fit(design[rows, ], prepared$treatment[rows],
                               family = stats::binomial())$coefficients
      direct[is.na(direct)] <- 0
      expect_equal(new$original_logit, as.numeric(design %*% direct))
      expect_equal(new$resample_logit, new$original_logit[rows])
    } else {
      refit <- stats::glm(rare_formula, data = resample, family = stats::binomial())
      expect_equal(new$original_logit,
                   as.numeric(stats::predict(refit, newdata = rare, type = "link")))
    }
  }
})

test_that("5d. bootstrap: a prediction failure, all levels present, uses the design rule", {
  # On real data predict() on a non-converged rank-deficient refit stops with
  # "NA/NaN/Inf in foreign function call"; here the failure is forced for chosen samples.
  prepared <- prepare_model_data(lalonde_formula, lalonde)
  n <- nrow(lalonde)
  failing <- c(2L, 5L)
  forced_error <- "NA/NaN/Inf in foreign function call (arg 1)"
  rows_of <- function(b) with_task_seed(derive_seed(7, 1000L + b), sample.int(n, n, replace = TRUE))
  failing_first_rows <- vapply(failing, function(b) rows_of(b)[1], integer(1))
  local_mocked_bindings(predict_bootstrap_refit = function(fit, data) {
    # The refit's own data identifies the sample: fail for the chosen samples only.
    if (nrow(fit$data) == n && any(vapply(failing_first_rows, function(r) {
      identical(fit$data$age[1], lalonde$age[r]) && identical(fit$data$re74[1], lalonde$re74[r])
    }, logical(1)))) stop(forced_error)
    as.numeric(stats::predict(fit, newdata = data, type = "link"))
  })
  ours <- match_each_model(model_set, lalonde, replace = TRUE, bootstrap_samples = 6, seed = 7)
  boot <- ours$bootstrap
  expect_equal(nrow(boot$skipped), 0)
  expect_equal(boot$prediction_fallback_samples, failing)
  expect_length(boot$design_rule_samples, 0)
  expect_false(anyNA(boot$rank))
  design <- cbind(1, prepared$design)
  for (b in seq_len(6)) {
    rows <- rows_of(b)
    resample <- lalonde[rows, , drop = FALSE]
    rownames(resample) <- NULL
    new <- bootstrap_sample_logits(lalonde_formula, lalonde, resample, rows, prepared)
    if (b %in% failing) {
      expect_identical(new$path, "design_fallback_predict_error")
      expect_false(new$design_rule)
      direct <- stats::glm.fit(design[rows, ], prepared$treatment[rows],
                               family = stats::binomial())$coefficients
      direct[is.na(direct)] <- 0
      expect_equal(new$original_logit, as.numeric(design %*% direct))
      expect_equal(new$resample_logit, new$original_logit[rows])
      expect_equal(rank(new$original_logit, ties.method = "average") / n, boot$rank[, b],
                   ignore_attr = TRUE)
    } else {
      expect_identical(new$path, "formula")
      refit <- stats::glm(lalonde_formula, data = resample, family = stats::binomial())
      expect_equal(new$original_logit,
                   as.numeric(stats::predict(refit, newdata = lalonde, type = "link")))
    }
  }
})

test_that("6. cores = 1 and cores = 4 give identical results", {
  # R CMD check --as-cran limits a package to 2 cores (_R_CHECK_LIMIT_CORES_); use 4 otherwise.
  parallel_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2 else 4
  for (m_order in c("largest", "random")) {
    serial <- match_each_model(model_set, lalonde, replace = TRUE, m_order = m_order,
                               oracle_logit = lalonde$age, bootstrap_samples = 4, cores = 1,
                               seed = 5)
    parallel <- match_each_model(model_set, lalonde, replace = TRUE, m_order = m_order,
                                 oracle_logit = lalonde$age, bootstrap_samples = 4,
                                 cores = parallel_cores, seed = 5)
    expect_equal(parallel$settings$cores, parallel_cores)
    expect_identical(serial$matches, parallel$matches)
    expect_identical(serial$weights, parallel$weights)
    expect_identical(serial$treated_kept, parallel$treated_kept)
    expect_identical(serial$oracle, parallel$oracle)
    expect_identical(serial$bootstrap, parallel$bootstrap)
  }
})

test_that("7. errors: wrong class, wrong row count, bad order, ratio other than 1", {
  expect_error(match_each_model(list(), lalonde), "propensity_model_set")
  expect_error(match_each_model(model_set, lalonde[-1, ]), "rows")
  shuffled <- lalonde[c(2:n_lalonde, 1), ]
  expect_error(match_each_model(model_set, shuffled), "treatment")
  expect_error(match_each_model(model_set, lalonde, m_order = "biggest"), "m_order")
  expect_error(match_each_model(model_set, lalonde, ratio = 2), "not supported in this version")
  expect_error(match_each_model(model_set, lalonde, replace = NA), "replace")
  expect_error(match_each_model(model_set, lalonde, caliper = -1), "caliper")
  expect_error(match_each_model(model_set, lalonde, oracle_logit = 1:3), "oracle_logit")
  expect_error(match_each_model(model_set, lalonde, bootstrap_samples = -1), "bootstrap_samples")
})

test_that("7b. common_random_order gives every model the same random order", {
  # Two models with identical scores: with one shared random order they get identical matches;
  # with an order per model (the default) the greedy matches without replacement differ.
  twin_set <- model_set
  twin_set$logit[, 2] <- twin_set$logit[, 1]
  common <- match_each_model(twin_set, lalonde, replace = FALSE, m_order = "random",
                             common_random_order = TRUE, seed = 3)
  expect_identical(common$matches[[1]]$control_id, common$matches[[2]]$control_id)
  expect_true(common$settings$common_random_order)
  separate <- match_each_model(twin_set, lalonde, replace = FALSE, m_order = "random", seed = 3)
  expect_false(identical(separate$matches[[1]]$control_id, separate$matches[[2]]$control_id))
  expect_false(separate$settings$common_random_order)
  expect_error(match_each_model(model_set, lalonde, common_random_order = TRUE),
               "needs `m_order = \"random\"`")
  expect_error(match_each_model(model_set, lalonde, m_order = "random",
                                common_random_order = NA), "TRUE or FALSE")
})

test_that("8. the stored matchit call re-evaluates to the same pairs and weights", {
  for (replace in c(FALSE, TRUE)) {
    ours <- match_each_model(model_set, lalonde, replace = replace, caliper = 0.5)
    for (id in ours$model_ids) {
      expect_match(ours$calls[[id]], "^MatchIt::matchit\\(")
      expect_match(ours$calls[[id]], "normalize = FALSE")
      expect_match(ours$calls[[id]], "std.caliper = TRUE")
      reference <- eval(parse(text = ours$calls[[id]]))
      matched <- !is.na(reference$match.matrix[, 1])
      expect_identical(ours$matches[[id]]$treated_id, rownames(reference$match.matrix)[matched])
      expect_identical(ours$matches[[id]]$control_id,
                       unname(reference$match.matrix[matched, 1]))
      expect_equal(unname(ours$weights[, id]), unname(reference$weights))
    }
  }
})

test_that("the object stores no data, and print lists every model with distinct controls", {
  ours <- match_each_model(model_set, lalonde, replace = TRUE)
  expect_s3_class(ours, "match_set")
  expect_false("data" %in% names(ours))
  expect_identical(ours$row_ids, rownames(lalonde))
  expect_identical(ours$treatment, as.integer(lalonde$treat))
  printed <- capture.output(print(ours))
  expect_true(any(grepl("distinct_controls", printed)))
  for (id in ours$model_ids) expect_true(any(grepl(id, printed, fixed = TRUE)))
  expect_true(any(grepl("with replacement", printed)))
})

test_that("9. timing: 7 models, both arms, B = 20, n = 2000 in under 10 seconds", {
  skip_unless_timing_tests()
  big <- simulate_example_data(n = 2000, scenario = "S2", seed = 1)
  big_set <- fit_propensity_models(treat ~ x1 + x2 + x3 + x4 + x5 + x6, big, n_trees = 200,
                                   algorithms = faster_set)
  expect_equal(ncol(big_set$logit), 7)
  timing <- system.time({
    with_replacement <- match_each_model(big_set, big, replace = TRUE, bootstrap_samples = 20,
                                         oracle_logit = big$true_logit)
    without_replacement <- match_each_model(big_set, big, replace = FALSE,
                                            bootstrap_samples = 20,
                                            oracle_logit = big$true_logit)
  })
  expect_elapsed_under(timing, 10, "match_each_model, 7 models, both arms, B = 20, n = 2000")
  expect_equal(dim(with_replacement$bootstrap$rank), c(2000, 20))
  expect_equal(length(without_replacement$matches), 7)
})
