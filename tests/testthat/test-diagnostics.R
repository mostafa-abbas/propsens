data("lalonde", package = "MatchIt")
lalonde_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75
lalonde_covariates <- lalonde[, c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
n_lalonde <- nrow(lalonde)
treated_lalonde <- lalonde$treat == 1
seeded_set <- fit_propensity_models(lalonde_formula, lalonde, seeds = 1:2, n_trees = 50,
                                    algorithms = faster_set)
logistic_set <- fit_propensity_models(lalonde_formula, lalonde,
                                      algorithms = "logistic_regression", cross_fit = FALSE)

# Three identical logistic columns under the three algorithm names (test 2).
identical_logit <- matrix(rep(logistic_set$logit[, 1], 3), ncol = 3)
identical_set <- fake_model_set(identical_logit,
                                c("logistic_regression", "random_forest", "gradient_boosting"),
                                lalonde$treat, formula = lalonde_formula)

test_that("1. invariance: an increasing transform of the logits leaves the rank measures", {
  before <- rank_disagreement(seeded_set)
  transformed <- seeded_set
  transformed$logit <- transformed$logit^3 + 2 * transformed$logit
  transformed$rank <- apply(transformed$logit, 2,
                            function(column) rank(column, ties.method = "average") / n_lalonde)
  after <- rank_disagreement(transformed)
  expect_equal(after$patients$rank_disagreement, before$patients$rank_disagreement)
  expect_equal(after$patients$rank_disagreement_between,
               before$patients$rank_disagreement_between)
  expect_equal(after$patients$rank_disagreement_within,
               before$patients$rank_disagreement_within)
  expect_false(isTRUE(all.equal(after$patients$logit_disagreement,
                                before$patients$logit_disagreement)))
})

test_that("2. identical models: every across-model diagnostic is at its no-disagreement value", {
  disagreement <- rank_disagreement(identical_set)
  # Zero up to floating point: on some platforms (macOS arm64) the SD of identical values is
  # about 1e-17 rather than exactly 0.
  expect_lt(max(abs(disagreement$patients$rank_disagreement)), 1e-12)
  expect_lt(max(abs(disagreement$patients$rank_disagreement_between)), 1e-12)
  expect_true(all(is.na(disagreement$patients$rank_disagreement_within)))
  expect_true(is.na(disagreement$summary$median_between_within_ratio_treated))
  for (replace in c(TRUE, FALSE)) {
    matched <- match_each_model(identical_set, lalonde, replace = replace)
    spread <- control_spread(matched, lalonde_covariates)
    expect_true(all(spread$patients$control_spread[treated_lalonde] == 0))
    expect_true(all(spread$patients$distinct_controls[treated_lalonde] == 1))
    balance <- balance_spread(matched, lalonde$treat, lalonde_covariates)
    expect_true(all(balance$spread$range == 0))
    expect_true(all(balance$spread$sd == 0))
    recurrence <- pair_recurrence(matched)
    expect_true(all(recurrence$patients$pair_recurrence_models[treated_lalonde] == 1))
    expect_equal(recurrence$summary$floor_models, 1 / 3)
  }
})

test_that("3. hand-computed 3-patient, 3-model example for d, d_between and d_within", {
  rank <- rbind(c(1, 2, 3), c(2, 1, 1), c(3, 3, 2)) / 3
  set <- fake_model_set(matrix(rnorm(9), 3),
                        c("logistic_regression", "random_forest", "random_forest"),
                        treatment = c(1, 0, 1), rank = rank)
  result <- rank_disagreement(set)
  expect_equal(result$patients$rank_disagreement, c(1 / 3, sqrt(3) / 9, sqrt(3) / 9))
  # Algorithm means: logistic = column 1, random forest = mean of columns 2 and 3.
  expect_equal(result$patients$rank_disagreement_between,
               c(1 / 2, 1 / 3, 1 / 6) / sqrt(2))
  # Only the random forest has two seeds: sqrt(variance of two values) = |difference| / sqrt(2).
  expect_equal(result$patients$rank_disagreement_within, c(1 / 3, 0, 1 / 3) / sqrt(2))
  expect_equal(result$summary$mean_d_treated, mean(c(1 / 3, sqrt(3) / 9)))
  expect_equal(result$summary$mean_d_all, mean(c(1 / 3, sqrt(3) / 9, sqrt(3) / 9)))
  expect_equal(result$summary$median_between_within_ratio_treated,
               stats::median(c((1 / 2) / (1 / 3), (1 / 6) / (1 / 3))))
  expect_equal(result$settings$algorithms_with_seeds, "random_forest")
})

test_that("4. control_spread: identical controls, hand-computed Mahalanobis, input errors", {
  covariates <- data.frame(x1 = c(0, 1, 2, 2, 5, 3, 1, 4), x2 = c(1, 0, 3, 3, 1, 4, 2, 0))
  treatment <- c(1, 1, 0, 0, 0, 0, 0, 0)
  pairs <- list(logistic_regression_1 = cbind(c(1, 2), c(3, 5)),
                random_forest_1 = cbind(c(1, 2), c(4, 6)))
  set <- fake_match_set(treatment, pairs, algorithms = c("logistic_regression", "random_forest"))
  spread <- control_spread(set, covariates)
  # Rows 3 and 4 have identical covariates: spread 0 although there are two distinct controls.
  expect_equal(spread$patients$distinct_controls[1], 2)
  expect_equal(spread$patients$control_spread[1], 0)
  covariance <- stats::cov(covariates)
  expected <- sqrt(stats::mahalanobis(unlist(covariates[5, ]), unlist(covariates[6, ]),
                                      covariance))
  expect_equal(spread$patients$control_spread[2], unname(expected))
  logistic_distance <- mean(sqrt(c(
    stats::mahalanobis(unlist(covariates[1, ]), unlist(covariates[3, ]), covariance),
    stats::mahalanobis(unlist(covariates[2, ]), unlist(covariates[5, ]), covariance)
  )))
  expect_equal(spread$summary$scale_distance, logistic_distance)
  expect_equal(spread$summary$scale_model, "logistic_regression_1")
  expect_equal(spread$patients$control_spread_scaled[2], unname(expected) / logistic_distance)
  expect_true(all(is.na(spread$patients$control_spread[treatment == 0])))
  duplicated_column <- cbind(covariates, x1_copy = covariates$x1)
  expect_message(with_copy <- control_spread(set, duplicated_column), "x1_copy")
  expect_equal(with_copy$patients$control_spread, spread$patients$control_spread)
  expect_message(control_spread(set, cbind(covariates, constant = 1)), "constant")
  expect_error(control_spread(set, covariates[-1, ]), "rows")
  expect_error(control_spread(set, cbind(covariates, label = letters[1:8])), "Character")
  expect_error(control_spread(list(), covariates), "match_set")
})

test_that("4b. control_spread: one shared control gives 0; three controls average 3 distances", {
  covariates <- data.frame(x1 = c(0, 1, 2, 7, 5, 3, 1, 4), x2 = c(1, 0, 3, 3, 1, 4, 2, 0))
  treatment <- c(1, 1, 0, 0, 0, 0, 0, 0)
  pairs <- list(a = cbind(c(1, 2), c(3, 3)), b = cbind(c(1, 2), c(3, 4)),
                c = cbind(c(1, 2), c(3, 5)))
  spread <- control_spread(fake_match_set(treatment, pairs), covariates)
  expect_equal(spread$patients$control_spread[1], 0)
  distances <- as.matrix(stats::dist(as.matrix(covariates) %*%
                                       solve(chol(stats::cov(covariates)))))
  expect_equal(spread$patients$control_spread[2],
               mean(c(distances[3, 4], distances[3, 5], distances[4, 5])))
  expect_equal(spread$summary$scale_model, "a")
})

test_that("5. balance_spread agrees with cobalt applied directly to the matchit object", {
  for (replace in c(TRUE, FALSE)) {
    matched <- match_each_model(logistic_set, lalonde, replace = replace)
    balance <- balance_spread(matched, lalonde$treat, lalonde_covariates)
    reference <- eval(parse(text = matched$calls[[1]]),
                      list(logistic_set = logistic_set, lalonde = lalonde))
    main <- cobalt::bal.tab(reference, binary = "std", un = FALSE)$Balance
    main <- main[main$Type != "Distance", ]
    expect_equal(balance$models$smd_max_main, max(abs(main$Diff.Adj)))
    extended <- cobalt::bal.tab(reference, binary = "std", int = TRUE, poly = 2,
                                un = FALSE)$Balance
    extended <- extended[extended$Type != "Distance", ]
    expect_equal(balance$models$smd_max_extended, max(abs(extended$Diff.Adj)))
    energy <- cobalt::bal.compute(cobalt::bal.init(lalonde_covariates, treat = lalonde$treat,
                                                   stat = "energy.dist", estimand = "ATT"),
                                  weights = reference$weights)
    ks <- cobalt::bal.compute(cobalt::bal.init(lalonde_covariates, treat = lalonde$treat,
                                               stat = "ks.max", estimand = "ATT"),
                              weights = reference$weights)
    expect_equal(balance$models$energy_distance, energy)
    expect_equal(balance$models$ks_max, ks)
    # bal.compute's smd.max standardizes binary covariates too, so it agrees.
    smd <- cobalt::bal.compute(cobalt::bal.init(lalonde_covariates, treat = lalonde$treat,
                                                stat = "smd.max", estimand = "ATT"),
                               weights = reference$weights)
    expect_equal(balance$models$smd_max_main, smd)
  }
})

test_that("5b. balance_spread: spread table, in_use column and input errors", {
  matched <- match_each_model(seeded_set, lalonde, replace = TRUE)
  balance <- balance_spread(matched, lalonde$treat, lalonde_covariates)
  expect_equal(balance$models$model_id, matched$model_ids)
  expect_equal(balance$models$in_use, matched$models$passed_logloss_gate)
  used <- balance$models[balance$models$in_use, ]
  for (statistic in balance$spread$statistic) {
    row <- balance$spread[balance$spread$statistic == statistic, ]
    expect_equal(row$range, max(used[[statistic]]) - min(used[[statistic]]))
    expect_equal(row$sd, stats::sd(used[[statistic]]))
    expect_equal(row$worst_value, max(used[[statistic]]))
    expect_equal(row$worst_model, used$model_id[which.max(used[[statistic]])])
  }
  expect_true(all(balance$models$smd_max_extended >= balance$models$smd_max_main))
  expect_error(balance_spread(matched, rev(lalonde$treat), lalonde_covariates), "treatment")
  expect_error(balance_spread(matched, lalonde$treat, lalonde_covariates[-1, ]), "rows")
  expect_error(balance_spread(matched, lalonde$treat, lalonde_covariates, polynomial = 0),
               "polynomial")
})

test_that("6. pair_recurrence: all-different controls give 1 / M", {
  treatment <- c(1, 1, 0, 0, 0, 0)
  pairs <- list(a = cbind(c(1, 2), c(3, 3)), b = cbind(c(1, 2), c(4, 3)),
                c = cbind(c(1, 2), c(5, 6)))
  recurrence <- pair_recurrence(fake_match_set(treatment, pairs))
  expect_equal(recurrence$patients$pair_recurrence_models[1], 1 / 3)
  expect_equal(recurrence$patients$pair_recurrence_models[2], 2 / 3)
  expect_equal(recurrence$patients$modal_control_models[2], "3")
  # Patient 1 has a three-way tie: the first control in the data is the modal one.
  expect_equal(recurrence$patients$modal_control_models[1], "3")
  expect_equal(recurrence$summary$floor_models, 1 / 3)
  expect_equal(recurrence$summary$bootstrap_samples_used, 0)
  expect_true(is.na(recurrence$summary$mean_bootstrap_treated))
})

test_that("6b. pair_recurrence bootstrap: corrected 1, uncorrected about 1 - 1/e", {
  n <- 50
  n_samples <- 200
  treatment <- c(1, rep(0, n - 1))
  present <- matrix(FALSE, n, n_samples)
  pairs <- matrix(NA_character_, n, n_samples)
  set.seed(11)
  for (b in seq_len(n_samples)) {
    rows <- sample.int(n, n, replace = TRUE)
    present[rows, b] <- TRUE
    if (present[1, b]) {
      # Control 2 is chosen whenever it is in the resample; otherwise one of the others in turn.
      others <- setdiff(which(present[, b]), 1:2)
      pairs[1, b] <- if (present[2, b]) "2" else as.character(others[(b %% length(others)) + 1])
    }
  }
  bootstrap <- list(rank = matrix(0.5, n, n_samples), pairs = pairs, present = present,
                    skipped = data.frame(sample = integer(0), reason = character(0)))
  set <- fake_match_set(treatment, list(a = cbind(1, 2), b = cbind(1, 2)),
                        bootstrap = bootstrap)
  recurrence <- pair_recurrence(set)
  expect_equal(recurrence$patients$modal_control_bootstrap[1], "2")
  expect_equal(recurrence$patients$pair_recurrence_bootstrap[1], 1)
  uncorrected <- recurrence$patients$pair_recurrence_bootstrap_uncorrected[1]
  expect_equal(uncorrected, sum(present[1, ] & present[2, ]) / sum(present[1, ]))
  expect_lt(abs(uncorrected - (1 - exp(-1))), 0.1)
  expect_equal(recurrence$summary$expected_ceiling_bootstrap_uncorrected, 1 - exp(-1))
  # A skipped resample (all-NA column) is ignored.
  bootstrap$present[, 1] <- NA
  bootstrap$pairs[, 1] <- NA
  skipped <- pair_recurrence(fake_match_set(treatment, list(a = cbind(1, 2)),
                                            bootstrap = bootstrap))
  expect_equal(skipped$summary$bootstrap_samples_used, n_samples - 1)
  expect_equal(skipped$patients$pair_recurrence_bootstrap[1], 1)
})

test_that("7. gate rule: one model passes, all models are used and every function flags it", {
  strict_set <- fit_propensity_models(lalonde_formula, lalonde, seeds = 1, n_trees = 50,
                                      quality_tolerance = 0, algorithms = faster_set)
  expect_equal(sum(strict_set$models$passed_logloss_gate), 1)
  matched <- match_each_model(strict_set, lalonde, replace = TRUE)
  results <- list(rank_disagreement(strict_set), control_spread(matched, lalonde_covariates),
                  balance_spread(matched, lalonde$treat, lalonde_covariates),
                  pair_recurrence(matched))
  for (result in results) {
    expect_true(result$quality_gate_left_fewer_than_two_models)
    expect_identical(result$models_used, strict_set$models$model_id)
    expect_output(print(result), "fewer than two")
  }
  normal <- list(rank_disagreement(seeded_set),
                 pair_recurrence(match_each_model(seeded_set, lalonde)))
  for (result in normal) {
    expect_false(result$quality_gate_left_fewer_than_two_models)
    expect_identical(result$models_used,
                     seeded_set$models$model_id[seeded_set$models$passed_logloss_gate])
  }
  ungated <- rank_disagreement(seeded_set, use_gate = FALSE)
  expect_identical(ungated$models_used, seeded_set$models$model_id)
  expect_error(rank_disagreement(seeded_set, use_gate = NA), "use_gate")
  expect_error(rank_disagreement(list()), "propensity_model_set")
})

test_that("8. the oracle entry and the bootstrap block never enter across-model statistics", {
  plain <- match_each_model(seeded_set, lalonde, replace = TRUE)
  wild_oracle <- -seeded_set$logit[, 1]
  extra <- match_each_model(seeded_set, lalonde, replace = TRUE, oracle_logit = wild_oracle,
                            bootstrap_samples = 5)
  expect_false(identical(extra$oracle$matches, plain$matches[[1]]))
  expect_identical(control_spread(extra, lalonde_covariates)$patients,
                   control_spread(plain, lalonde_covariates)$patients)
  expect_identical(balance_spread(extra, lalonde$treat, lalonde_covariates)$models,
                   balance_spread(plain, lalonde$treat, lalonde_covariates)$models)
  recurrence_extra <- pair_recurrence(extra)
  recurrence_plain <- pair_recurrence(plain)
  expect_identical(recurrence_extra$patients$pair_recurrence_models,
                   recurrence_plain$patients$pair_recurrence_models)
  expect_false("oracle" %in% recurrence_extra$models_used)
  expect_equal(recurrence_extra$summary$bootstrap_samples_used, 5)
})

test_that("9. timing: the four diagnostics on n = 2000 and 7 models in under 5 seconds", {
  skip_unless_timing_tests()
  big <- simulate_example_data(n = 2000, scenario = "S2", seed = 1)
  big_set <- fit_propensity_models(treat ~ x1 + x2 + x3 + x4 + x5 + x6, big, n_trees = 200,
                                   algorithms = faster_set)
  big_matches <- match_each_model(big_set, big, replace = TRUE, bootstrap_samples = 20)
  covariates <- big[, paste0("x", 1:6)]
  expect_equal(ncol(big_set$logit), 7)
  timing <- system.time({
    rank_disagreement(big_set)
    control_spread(big_matches, covariates)
    balance_spread(big_matches, big$treat, covariates)
    pair_recurrence(big_matches)
  })
  expect_elapsed_under(timing, 5, "four diagnostics, 7 models, n = 2000")
})

test_that("the extended max |SMD| over covariate block pairs equals the full table", {
  set.seed(84)
  n <- 300
  covariates <- as.data.frame(matrix(stats::rnorm(n * 18), n, 18))
  covariates$b1 <- stats::rbinom(n, 1, 0.3)
  covariates$b2 <- stats::rbinom(n, 1, 0.6)
  covariates$f1 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  covariates$f2 <- factor(sample(c("u", "v", "w", "z"), n, replace = TRUE))
  covariates$constant_product <- 0   # all-zero: its terms are dropped by cobalt in every table
  treatment <- stats::rbinom(n, 1, 0.3)
  weights <- list(m1 = ifelse(treatment == 1, 1, stats::runif(n)),
                  m2 = ifelse(treatment == 1, 1, stats::rbinom(n, 2, 0.4)))
  statistics <- suppressWarnings(balance_statistics(covariates, treatment, weights))
  full <- suppressWarnings(vapply(names(weights), function(id) {
    table <- cobalt::bal.tab(covariates, treat = treatment, weights = weights, method = "matching",
                             estimand = "ATT", s.d.denom = "treated", binary = "std",
                             continuous = "std", int = TRUE, poly = 2, un = FALSE,
                             quick = TRUE)$Balance
    max(abs(table[[paste0("Diff.", id)]]), na.rm = TRUE)
  }, numeric(1)))
  expect_gt(ncol(covariates), 20)   # the block path is used
  expect_identical(statistics$smd_max_extended, unname(full))
})
