# logistic regression variants in fit_propensity_models().
variants_formula <- treat ~ x1 + x2 + x3 + x4 + x5 + x6
variants_data <- simulate_example_data(n = 300, scenario = "S2", seed = 6)
variants_design <- stats::model.matrix(variants_formula, variants_data)[, -1]

# Independent re-implementations for the leakage test.
hand_splines <- function(x_train, treatment, x_test) {
  splined <- apply(x_train, 2, function(v) {
    knots <- stats::quantile(v, c(1, 2) / 3, names = FALSE)
    length(unique(v)) >= 5 && all(diff(c(min(v), knots, max(v))) > 0)
  })
  expand <- function(x) {
    do.call(cbind, lapply(seq_len(ncol(x)), function(j) {
      if (splined[j]) {
        stats::predict(splines::ns(x_train[, j], df = 3), x[, j])
      } else {
        x[, j, drop = FALSE]
      }
    }))
  }
  fit <- stats::glm.fit(cbind(1, expand(x_train)), treatment, family = stats::binomial())
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  as.numeric(stats::plogis(cbind(1, expand(x_test)) %*% beta))
}
hand_lasso <- function(x_train, treatment, x_test) {
  folds <- sample(rep_len(1:5, nrow(x_train)))
  fit <- glmnet::cv.glmnet(x_train, treatment, family = "binomial", foldid = folds, alpha = 1)
  as.numeric(stats::predict(fit, newx = x_test, s = "lambda.min", type = "response"))
}

test_that("the logistic variants are not in the default set", {
  defaults <- eval(formals(fit_propensity_models)$algorithms)
  expect_false(any(c("logistic_splines", "lasso", "lasso_interactions") %in% defaults))
})

test_that("no leakage, a hand-written cross-fitting loop reproduces each variant", {
  skip_if_not_installed("glmnet")
  interactions <- two_way_interaction_design(variants_design)
  for (algorithm in c("logistic_splines", "lasso", "lasso_interactions")) {
    set <- fit_propensity_models(variants_formula, variants_data, algorithms = algorithm,
                                 seeds = 1, folds = 5, seed = 1)
    x <- if (algorithm == "lasso_interactions") interactions else variants_design
    hand <- numeric(nrow(x))
    hand_loop <- function() {
      for (k in 1:5) {
        train <- set$folds != k
        hand[!train] <<- if (algorithm == "logistic_splines") {
          hand_splines(x[train, ], variants_data$treat[train], x[!train, ])
        } else {
          hand_lasso(x[train, ], variants_data$treat[train], x[!train, ])
        }
      }
    }
    withr::with_seed(set$models$task_seed[1], hand_loop(), .rng_kind = "Mersenne-Twister",
                     .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    hand <- pmin(pmax(hand, 0.001), 0.999)
    expect_equal(unname(set$probability[, 1]), hand, tolerance = 1e-10, info = algorithm)
  }
})

test_that("same seed identical, seeds differ for the lasso, cores 1 vs 4 identical", {
  skip_if_not_installed("glmnet")
  parallel_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2 else 4
  algorithms <- c("logistic_regression", "logistic_splines", "lasso", "lasso_interactions")
  run <- function(cores) {
    fit_propensity_models(variants_formula, variants_data, algorithms = algorithms,
                          seeds = 1:2, seed = 8, cores = cores)
  }
  first <- run(1)
  expect_identical(first$probability, run(1)$probability)
  expect_identical(first$probability, run(parallel_cores)$probability)
  expect_identical(colnames(first$probability),
                   c("logistic_regression_1", "logistic_splines_1", "lasso_1", "lasso_2",
                     "lasso_interactions_1", "lasso_interactions_2"))
  expect_false(isTRUE(all.equal(first$probability[, "lasso_interactions_1"],
                                first$probability[, "lasso_interactions_2"])))
})

test_that("clipping keeps the variants inside [clip, 1 - clip]", {
  skip_if_not_installed("glmnet")
  separable <- data.frame(treat = rep(0:1, each = 60), x1 = c(rnorm(60, -5), rnorm(60, 5)),
                          x2 = rnorm(120))
  set <- suppressWarnings(fit_propensity_models(
    treat ~ x1 + x2, separable, seeds = 1, clip = 0.01, seed = 2,
    algorithms = c("logistic_splines", "lasso", "lasso_interactions")
  ))
  expect_true(all(is.finite(set$logit)))
  expect_true(all(set$probability >= 0.01 & set$probability <= 0.99))
})

test_that("splines enter only for columns with 5+ values and distinct interior knots", {
  concentrated <- c(rep(0, 30), 1:10)
  expect_false(spline_eligible(concentrated, 3L, 5L))
  expect_true(spline_eligible(seq_len(40), 3L, 5L))
  x <- cbind(continuous = seq_len(40) / 7, binary = rep(0:1, 20))
  treatment <- rep(c(0, 1, 1, 0), 10)
  with_splines <- fit_logistic_splines(x, treatment, x,
                                       list(df = 3L, minimum_distinct_values = 5L))
  linear_only <- fit_logistic_splines(x, treatment, x,
                                      list(df = 3L, minimum_distinct_values = 100L))
  expect_false(isTRUE(all.equal(with_splines, linear_only)))
  expect_equal(linear_only, as.numeric(stats::fitted(stats::glm(treatment ~ x,
                                                                family = binomial()))),
               tolerance = 1e-8)
})
