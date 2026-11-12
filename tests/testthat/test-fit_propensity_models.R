example_formula <- treat ~ x1 + x2 + x3 + x4 + x5 + x6
example_data <- simulate_example_data(n = 600, scenario = "S2", seed = 1)
quick_set <- function(data = example_data, algorithms = faster_set, ...) {
  # nolint next: object_usage_linter. The package function is attached by testthat.
  fit_propensity_models(example_formula, data, n_trees = 100, algorithms = algorithms, ...)
}

test_that("no leakage: a hand-written cross-fitting loop reproduces the logistic column", {
  set <- quick_set(algorithms = "logistic_regression", folds = 5, seed = 1)
  hand <- numeric(nrow(example_data))
  for (k in 1:5) {
    train <- set$folds != k
    fit <- glm(example_formula, data = example_data[train, ], family = binomial())
    hand[!train] <- predict(fit, newdata = example_data[!train, ], type = "response")
  }
  hand <- pmin(pmax(hand, 0.001), 0.999)
  expect_equal(unname(set$logit[, "logistic_regression_1"]), log(hand / (1 - hand)),
               tolerance = 1e-12)
  expect_equal(sort(unique(set$folds)), 1:5)
  for (group in 0:1) {
    counts <- table(set$folds[example_data$treat == group])
    expect_lte(diff(range(counts)), 1)
  }
})

test_that("reproducible: same seed twice, and cores = 1 versus cores = 4, identical", {
  # R CMD check --as-cran limits a package to 2 cores (_R_CHECK_LIMIT_CORES_); use 4 otherwise.
  parallel_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2 else 4
  serial_1 <- quick_set(seeds = 1:2, seed = 3, cores = 1)
  serial_2 <- quick_set(seeds = 1:2, seed = 3, cores = 1)
  parallel_4 <- quick_set(seeds = 1:2, seed = 3, cores = parallel_cores)
  expect_equal(parallel_4$settings$cores, parallel_cores)
  expect_identical(serial_1$logit, serial_2$logit)
  expect_identical(serial_1$rank, serial_2$rank)
  expect_identical(serial_1$probability, serial_2$probability)
  expect_identical(serial_1$logit, parallel_4$logit)
  expect_identical(serial_1$rank, parallel_4$rank)
  expect_identical(serial_1$models, parallel_4$models)
})

test_that("different seeds give different tree columns; logistic regression has one column", {
  set <- quick_set(seeds = 1:3, seed = 1)
  expect_identical(colnames(set$logit),
                   c("logistic_regression_1", paste0("random_forest_", 1:3),
                     paste0("gradient_boosting_", 1:3)))
  expect_equal(sum(set$models$algorithm == "logistic_regression"), 1)
  expect_false(isTRUE(all.equal(set$logit[, "random_forest_1"], set$logit[, "random_forest_2"])))
  expect_false(isTRUE(all.equal(set$logit[, "gradient_boosting_1"],
                                set$logit[, "gradient_boosting_2"])))
  expect_equal(length(unique(set$models$task_seed)), nrow(set$models))
})

test_that("clipping: no infinite logit on a perfectly separable dataset", {
  separable <- data.frame(treat = rep(0:1, each = 60), x1 = c(rnorm(60, -5), rnorm(60, 5)),
                          x2 = rnorm(120))
  set <- suppressWarnings(fit_propensity_models(treat ~ x1 + x2, separable, seeds = 1,
                                                n_trees = 50, clip = 0.01, seed = 2,
                                                algorithms = faster_set))
  expect_true(all(is.finite(set$logit)))
  expect_true(all(set$probability >= 0.01 & set$probability <= 0.99))
  expect_true(any(set$probability == 0.01) || any(set$probability == 0.99))
})

test_that("ranks lie in [1/n, 1] and are monotone in the logits within each model", {
  set <- quick_set(seeds = 1, seed = 1)
  n <- nrow(example_data)
  expect_true(all(set$rank >= 1 / n & set$rank <= 1))
  for (id in colnames(set$logit)) {
    expect_equal(set$rank[, id], rank(set$logit[, id], ties.method = "average") / n)
    ordered <- order(set$logit[, id])
    expect_true(all(diff(set$rank[ordered, id]) >= 0))
  }
})

test_that("log-loss, Brier and c-statistic agree with hand-computed values", {
  probability <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  treatment <- c(0, 0, 1, 0, 0, 1, 1, 0, 1, 1)
  quality <- propsens:::model_quality(probability, treatment)
  expect_equal(unname(quality["logloss"]),
               -mean(treatment * log(probability) + (1 - treatment) * log(1 - probability)))
  expect_equal(unname(quality["brier"]), mean((probability - treatment)^2))
  pairs <- expand.grid(treated = probability[treatment == 1], control = probability[treatment == 0])
  expect_equal(unname(quality["c_statistic"]),
               mean((pairs$treated > pairs$control) + 0.5 * (pairs$treated == pairs$control)))
  expect_equal(propsens:::compute_c_statistic(rep(0.3, 10), treatment), 0.5)
  expect_equal(propsens:::compute_c_statistic(treatment, treatment), 1)
  tiny <- data.frame(treat = treatment, x1 = probability, x2 = c(1, 3, 2, 5, 4, 7, 6, 9, 8, 10))
  set <- fit_propensity_models(treat ~ x1 + x2, tiny, algorithms = "logistic_regression",
                               folds = 2, seed = 1)
  expected <- propsens:::model_quality(set$probability[, 1], treatment)
  expect_equal(set$models$logloss, unname(expected["logloss"]))
  expect_equal(set$models$brier, unname(expected["brier"]))
  expect_equal(set$models$c_statistic, unname(expected["c_statistic"]))
})

test_that("log-loss gate: tolerance 0 keeps only the best model, tolerance 1 keeps all", {
  strict <- quick_set(seeds = 1:2, seed = 1, quality_tolerance = 0)
  expect_equal(sum(strict$models$passed_logloss_gate), 1)
  expect_true(strict$models$passed_logloss_gate[which.min(strict$models$logloss)])
  loose <- quick_set(seeds = 1:2, seed = 1, quality_tolerance = 1)
  expect_true(all(loose$models$passed_logloss_gate))
  expect_identical(strict$logit, loose$logit)
})

test_that("n_trees: chosen once by cross-validation when NULL, used unchanged when supplied", {
  chosen <- fit_propensity_models(example_formula, example_data,
                                  algorithms = "gradient_boosting", seeds = 1:2, seed = 1)
  expect_true(is.integer(chosen$settings$n_trees))
  expect_gte(chosen$settings$n_trees, 1)
  expect_equal(chosen$settings$n_trees_source, "cross-validation")
  expect_true(all(chosen$models$n_trees == chosen$settings$n_trees))
  supplied <- fit_propensity_models(example_formula, example_data,
                                    algorithms = "gradient_boosting", seeds = 1:2, seed = 1,
                                    n_trees = 123)
  expect_equal(supplied$settings$n_trees, 123)
  expect_equal(supplied$settings$n_trees_source, "supplied")
  expect_true(all(supplied$models$n_trees == 123L))
  logistic_only <- quick_set(algorithms = "logistic_regression", seed = 1)
  expect_true(all(is.na(logistic_only$models$n_trees)))
})

test_that("errors are clear: treatment not 0/1, missing values, folds < 2", {
  bad_treatment <- example_data
  bad_treatment$treat <- bad_treatment$treat + 1
  expect_error(quick_set(bad_treatment), "coded 0/1")
  factor_treatment <- example_data
  factor_treatment$treat <- factor(factor_treatment$treat)
  expect_error(quick_set(factor_treatment), "coded 0/1")
  with_missing <- example_data
  with_missing$x2[5] <- NA
  expect_error(quick_set(with_missing), "[Ii]mpute first")
  expect_error(quick_set(folds = 1), "`folds`")
  expect_error(quick_set(seeds = integer(0)), "`seeds`")
  expect_error(fit_propensity_models(~ x1, example_data), "two-sided")
  expect_error(quick_set(clip = 0.5), "`clip`")
})

test_that("cross_fit = FALSE (as MatchIt) differs from cross-fitting for the tree models", {
  fitted <- quick_set(seeds = 1, seed = 1, cross_fit = TRUE)
  in_sample <- quick_set(seeds = 1, seed = 1, cross_fit = FALSE)
  expect_false(isTRUE(all.equal(fitted$logit[, "random_forest_1"],
                                in_sample$logit[, "random_forest_1"])))
  expect_false(isTRUE(all.equal(fitted$logit[, "gradient_boosting_1"],
                                in_sample$logit[, "gradient_boosting_1"])))
  expect_false(isTRUE(all.equal(fitted$logit[, "logistic_regression_1"],
                                in_sample$logit[, "logistic_regression_1"])))
  expect_false(in_sample$settings$cross_fit)
  expect_match(paste(capture.output(print(in_sample)), collapse = " "), "as MatchIt does it")
  rd_file <- testthat::test_path("..", "..", "man", "fit_propensity_models.Rd")
  if (file.exists(rd_file)) {
    rd_text <- paste(readLines(rd_file), collapse = " ")
  } else {
    rd_db <- tryCatch(tools::Rd_db("propsens"), error = function(e) NULL)
    skip_if(is.null(rd_db), "Rd database not available")
    rd_text <- paste(capture.output(tools::Rd2txt(rd_db[["fit_propensity_models.Rd"]])),
                     collapse = " ")
  }
  expect_match(rd_text, "MatchIt")
  expect_match(rd_text, "out-of-bag")
})

test_that("timing: n = 2,000, Scenario 2, default settings, single core, under 30 seconds", {
  skip_unless_timing_tests()
  big <- simulate_example_data(n = 2000, scenario = "S2", seed = 1)
  timing <- system.time(set <- fit_propensity_models(example_formula, big, cores = 1,
                                                     algorithms = faster_set))
  expect_elapsed_under(timing, 30, "fit_propensity_models, 7 models, n = 2000")
  expect_equal(dim(set$logit), c(2000, 7))
})

test_that("a factor level missing from a training fold does not stop cross-fitting (real data)", {
  set.seed(11)
  n <- 300
  frame <- data.frame(x1 = stats::rnorm(n),
                      x2 = factor(sample(c("common", "other"), n, replace = TRUE),
                                  levels = c("common", "other", "rare")))
  frame$x2[7] <- "rare"
  frame$treat <- stats::rbinom(n, 1, stats::plogis(-1 + frame$x1))
  # The rare level is in one row, so the fold holding it is predicted from training rows that
  # lack the level.
  fit <- suppressWarnings(fit_propensity_models(treat ~ x1 + x2, frame,
                                                algorithms = "logistic_regression", seeds = 1,
                                                seed = 3))
  expect_true(all(is.finite(fit$probability)))
  folds <- fit$folds
  train <- folds != folds[7]
  expect_true(has_unseen_levels(treat ~ x1 + x2, frame, train, !train))
  expect_false(has_unseen_levels(treat ~ x1, frame, train, !train))
  # Without unseen levels the formula path is used and matches glm() exactly.
  plain <- fit_propensity_models(treat ~ x1, frame, algorithms = "logistic_regression",
                                 seeds = 1, seed = 3)
  by_hand <- numeric(n)
  for (k in sort(unique(plain$folds))) {
    held <- plain$folds == k
    model <- stats::glm(treat ~ x1, data = frame[!held, ], family = stats::binomial())
    by_hand[held] <- stats::predict(model, newdata = frame[held, ], type = "response")
  }
  expect_equal(unname(plain$probability[, "logistic_regression_1"]),
               pmin(pmax(by_hand, 0.001), 0.999))
})

test_that("the design-matrix fallback predicts as predict() on a rank-deficient logistic fit", {
  # predict() on a rank-deficient fit that did not converge can stop (non-finite QR);
  # predict_held_out() then uses logistic_on_design(). Where predict() works, both agree.
  set.seed(12)
  n <- 400
  frame <- data.frame(x1 = stats::rnorm(n), x2 = factor(sample(letters[1:4], n, replace = TRUE)))
  frame$x3 <- 2 * frame$x1   # aliased: coefficient not estimable
  frame$treat <- stats::rbinom(n, 1, stats::plogis(-0.5 + frame$x1))
  formula <- treat ~ x1 + x2 + x3
  prepared <- prepare_model_data(formula, frame)
  train <- seq_len(n) <= 300
  test <- !train
  model <- stats::glm(formula, data = frame[train, ], family = stats::binomial())
  expected <- suppressWarnings(stats::predict(model, newdata = frame[test, ], type = "response"))
  expect_equal(logistic_on_design(prepared, train, test), unname(expected), tolerance = 1e-6)
})

test_that("logistic_log records the prediction path and glm's warnings per fold", {
  set.seed(11)
  n <- 300
  frame <- data.frame(x1 = stats::rnorm(n),
                      x2 = factor(sample(c("common", "other"), n, replace = TRUE),
                                  levels = c("common", "other", "rare")))
  frame$x2[7] <- "rare"
  frame$treat <- stats::rbinom(n, 1, stats::plogis(-1 + frame$x1))
  fit <- suppressWarnings(fit_propensity_models(treat ~ x1 + x2, frame,
                                                algorithms = "logistic_regression", seeds = 1,
                                                seed = 3))
  log <- fit$logistic_log
  expect_equal(names(log), c("model_id", "fold", "path", "glm_not_converged",
                             "glm_fitted_0_or_1"))
  expect_equal(nrow(log), fit$settings$folds)
  expect_equal(log$path[log$fold == fit$folds[7]], "design_unseen_levels")
  expect_true(all(log$path[log$fold != fit$folds[7]] == "formula"))
  # Complete separation: glm warns, and the warning is recorded.
  separated <- data.frame(x1 = c(stats::rnorm(150, -3), stats::rnorm(150, 3)),
                          treat = rep(0:1, each = 150))
  fit_separated <- suppressWarnings(fit_propensity_models(treat ~ x1, separated,
                                                          algorithms = "logistic_regression",
                                                          seeds = 1, seed = 3))
  expect_true(any(fit_separated$logistic_log$glm_fitted_0_or_1))
  # No logistic regression: no log.
  expect_null(suppressWarnings(fit_propensity_models(treat ~ x1, frame, algorithms = "bart",
                                                     seeds = 1, seed = 3))$logistic_log)
})
