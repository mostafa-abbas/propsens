# BART (dbarts::bart2) and the GRF probability forest in fit_propensity_models().
# Small settings keep the tests fast; the fitting code path is the one used with the defaults.
bart_grf_formula <- treat ~ x1 + x2 + x3 + x4 + x5 + x6
bart_grf_data <- simulate_example_data(n = 300, scenario = "S2", seed = 4)
small_settings <- function() {
  list(bart = list(n.trees = 20L, n.burn = 50L, n.samples = 50L, n.chains = 2L, n.threads = 1L),
       grf_probability = list(num.trees = 200L, honesty = TRUE, num.threads = 1L))
}
# The package's model fit with smaller BART/GRF settings (the same internal code path).
small_set <- function(algorithms, ...) {
  settings <- small_settings()
  testthat::local_mocked_bindings(
    fit_bart = function(x_train, treatment, x_test, ignored) {
      bart_fit(x_train, treatment, x_test, settings$bart)
    },
    fit_grf_probability = function(x_train, treatment, x_test, ignored) {
      grf_fit(x_train, treatment, x_test, settings$grf_probability)
    }
  )
  fit_propensity_models(bart_grf_formula, bart_grf_data, algorithms = algorithms, ...)
}
# Direct calls of the libraries, written independently of the package helpers.
bart_fit <- function(x_train, treatment, x_test, settings) {
  library_seed <- sample.int(.Machine$integer.max, 1L)
  fit <- dbarts::bart2(x_train, treatment, test = x_test, n.trees = settings$n.trees,
                       n.burn = settings$n.burn, n.samples = settings$n.samples,
                       n.chains = settings$n.chains, n.threads = 1L, seed = library_seed,
                       verbose = FALSE)
  as.numeric(stats::fitted(fit, type = "ev", sample = if (is.null(x_test)) "train" else "test"))
}
grf_fit <- function(x_train, treatment, x_test, settings) {
  library_seed <- sample.int(.Machine$integer.max, 1L)
  fit <- grf::probability_forest(x_train, factor(treatment, levels = c(0, 1)),
                                 num.trees = settings$num.trees, honesty = TRUE,
                                 num.threads = 1L, seed = library_seed)
  predictions <- if (is.null(x_test)) predict(fit)$predictions else
    predict(fit, newdata = x_test)$predictions
  as.numeric(predictions[, "1"])
}

test_that("the default set is the same as in propensity_sensitivity()", {
  expect_identical(eval(formals(fit_propensity_models)$algorithms),
                   c("logistic_regression", "bart", "grf_probability"))
  expect_identical(eval(formals(fit_propensity_models)$algorithms),
                   eval(formals(propensity_sensitivity)$algorithms))
})

test_that("no leakage, a hand-written cross-fitting loop reproduces BART and GRF", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("grf")
  settings <- small_settings()
  for (algorithm in c("bart", "grf_probability")) {
    set <- small_set(algorithm, seeds = 1, folds = 5, seed = 1)
    design <- stats::model.matrix(bart_grf_formula, bart_grf_data)[, -1]
    treatment <- bart_grf_data$treat
    spec_seed <- set$models$task_seed[1]
    hand <- numeric(nrow(design))
    hand_loop <- function() {
      for (k in 1:5) {
        train <- set$folds != k
        hand[!train] <<- if (algorithm == "bart") {
          bart_fit(design[train, ], treatment[train], design[!train, ], settings$bart)
        } else {
          grf_fit(design[train, ], treatment[train], design[!train, ], settings$grf_probability)
        }
      }
    }
    withr::with_seed(spec_seed, hand_loop(), .rng_kind = "Mersenne-Twister",
                     .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    hand <- pmin(pmax(hand, 0.001), 0.999)
    expect_equal(unname(set$probability[, 1]), hand, tolerance = 1e-12, info = algorithm)
  }
})

test_that("same seed twice identical, different seeds differ, cores 1 vs 4 identical", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("grf")
  parallel_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2 else 4
  algorithms <- c("logistic_regression", "bart", "grf_probability")
  first <- small_set(algorithms, seeds = 1:2, seed = 7, cores = 1)
  second <- small_set(algorithms, seeds = 1:2, seed = 7, cores = 1)
  parallel <- small_set(algorithms, seeds = 1:2, seed = 7, cores = parallel_cores)
  expect_identical(first$probability, second$probability)
  expect_identical(first$probability, parallel$probability)
  expect_false(isTRUE(all.equal(first$probability[, "bart_1"], first$probability[, "bart_2"])))
  expect_false(isTRUE(all.equal(first$probability[, "grf_probability_1"],
                                first$probability[, "grf_probability_2"])))
  expect_identical(colnames(first$probability),
                   c("logistic_regression_1", "bart_1", "bart_2", "grf_probability_1",
                     "grf_probability_2"))
})

test_that("clipping keeps BART and GRF probabilities inside [clip, 1 - clip]", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("grf")
  separable <- data.frame(treat = rep(0:1, each = 60), x1 = c(rnorm(60, -5), rnorm(60, 5)),
                          x2 = rnorm(120))
  settings <- small_settings()
  testthat::local_mocked_bindings(
    fit_bart = function(x_train, treatment, x_test, ignored) {
      bart_fit(x_train, treatment, x_test, settings$bart)
    },
    fit_grf_probability = function(x_train, treatment, x_test, ignored) {
      grf_fit(x_train, treatment, x_test, settings$grf_probability)
    }
  )
  set <- fit_propensity_models(treat ~ x1 + x2, separable, seeds = 1,
                               algorithms = c("bart", "grf_probability"), clip = 0.01, seed = 2)
  expect_true(all(is.finite(set$logit)))
  expect_true(all(set$probability >= 0.01 & set$probability <= 0.99))
})

test_that("cross_fit = FALSE runs for BART and GRF (in-sample / out-of-bag)", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("grf")
  set <- small_set(c("bart", "grf_probability"), seeds = 1, cross_fit = FALSE, seed = 3)
  expect_true(all(set$probability > 0 & set$probability < 1))
})

test_that("the package's BART and GRF helpers equal direct library calls", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("grf")
  settings <- small_settings()
  design <- stats::model.matrix(bart_grf_formula, bart_grf_data)[, -1]
  train <- seq_len(nrow(design)) <= 200
  for (pair in list(list(fit_bart, bart_fit, settings$bart),
                    list(fit_grf_probability, grf_fit, settings$grf_probability))) {
    for (x_test in list(design[!train, ], NULL)) {
      package_value <- withr::with_seed(5, pair[[1]](design[train, ], bart_grf_data$treat[train],
                                                     x_test, pair[[3]]))
      direct_value <- withr::with_seed(5, pair[[2]](design[train, ], bart_grf_data$treat[train],
                                                    x_test, pair[[3]]))
      expect_identical(package_value, direct_value)
    }
  }
})

test_that("propensity_sensitivity() defaults to logistic regression, BART and GRF", {
  expect_identical(eval(formals(propensity_sensitivity)$algorithms),
                   c("logistic_regression", "bart", "grf_probability"))
  expect_identical(eval(formals(propensity_sensitivity)$seeds), 1:3)
  settings <- small_settings()
  testthat::local_mocked_bindings(
    fit_bart = function(x_train, treatment, x_test, ignored) {
      bart_fit(x_train, treatment, x_test, settings$bart)
    },
    fit_grf_probability = function(x_train, treatment, x_test, ignored) {
      grf_fit(x_train, treatment, x_test, settings$grf_probability)
    }
  )
  result <- sensitivity_quietly(bart_grf_formula, bart_grf_data, bootstrap_samples = 2,
                                replace = FALSE, seed = 2)
  expect_identical(result$settings$algorithms, c("logistic_regression", "bart", "grf_probability"))
  expect_identical(result$models$model_id,
                   c("logistic_regression_1", "bart_1", "bart_2", "bart_3", "grf_probability_1",
                     "grf_probability_2", "grf_probability_3"))
})

test_that("the fixed BART and GRF settings of the default are as documented", {
  fitted_set <- small_set(c("logistic_regression", "bart", "grf_probability"), seeds = 1)
  expect_identical(fitted_set$settings$tree_settings$bart$n.trees, 75L)
  expect_identical(fitted_set$settings$tree_settings$bart$n.burn, 500L)
  expect_identical(fitted_set$settings$tree_settings$bart$n.samples, 500L)
  expect_identical(fitted_set$settings$tree_settings$bart$n.chains, 4L)
  expect_identical(fitted_set$settings$tree_settings$grf_probability$num.trees, 2000L)
  expect_true(fitted_set$settings$tree_settings$grf_probability$honesty)
})
