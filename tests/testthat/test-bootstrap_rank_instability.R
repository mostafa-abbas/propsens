# bootstrap_rank_instability(), as defined in the paper: the SD of a patient's percentile
# rank over the logistic bootstrap refits, averaged over treated patients.
data("lalonde", package = "MatchIt")
instability_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75

test_that("the value is the mean over treated of the SD of the rank; skipped samples drop out", {
  # Two usable samples and one skipped (a column of NA).
  ranks <- cbind(c(0.1, 0.2, 0.3, 0.4), c(0.3, 0.2, 0.5, 0.8), NA_real_)
  treatment <- c(1, 1, 0, 0)
  result <- instability_from_ranks(ranks, treatment, as.character(1:4), "test")
  expect_equal(result$patients$rank_instability, c(sqrt(0.02), 0, sqrt(0.02), sqrt(0.08)))
  expect_equal(result$summary$mean_treated, sqrt(0.02) / 2)
  expect_identical(result$summary$samples_used, 2L)
  expect_identical(result$summary$samples_skipped, 1L)
})

test_that("fewer than two usable samples give NA", {
  ranks <- cbind(c(0.1, 0.2), NA_real_)
  result <- instability_from_ranks(ranks, c(1, 0), c("a", "b"), "test")
  expect_true(is.na(result$summary$mean_treated))
})

test_that("a formula, a matchit object and stored ranks give the same value for the same seed", {
  model_set <- fit_propensity_models(instability_formula, lalonde,
                                     algorithms = "logistic_regression", seed = 3)
  match_set <- match_each_model(model_set, lalonde, replace = FALSE, bootstrap_samples = 5,
                                seed = 3)
  stored <- bootstrap_rank_instability(match_set)
  from_formula <- bootstrap_rank_instability(instability_formula, lalonde, bootstrap_samples = 5,
                                             seed = 3)
  expect_equal(from_formula$summary$mean_treated, stored$summary$mean_treated)
  expect_equal(from_formula$patients, stored$patients)
  ranks <- match_set$bootstrap$rank
  expect_equal(stored$summary$mean_treated,
               mean(apply(ranks[lalonde$treat == 1, ], 1, stats::sd)))
  matched <- MatchIt::matchit(instability_formula, data = lalonde, distance = "glm",
                              link = "linear.logit")
  from_matchit <- bootstrap_rank_instability(matched, bootstrap_samples = 5, seed = 3)
  expect_equal(from_matchit$summary$mean_treated, stored$summary$mean_treated)
})

test_that("print states the value and the number of refits", {
  result <- bootstrap_rank_instability(instability_formula, lalonde, bootstrap_samples = 3)
  printed <- utils::capture.output(print(result))
  expect_match(printed[1], "Bootstrap rank instability of the logistic model")
  expect_match(printed[1], format_number(result$summary$mean_treated), fixed = TRUE)
  expect_match(printed[2], "3 bootstrap refits")
})

test_that("errors are clear", {
  expect_error(bootstrap_rank_instability("treat"), "propensity_sensitivity object")
  expect_error(bootstrap_rank_instability(instability_formula), "pass the data")
  expect_error(bootstrap_rank_instability(instability_formula, lalonde, bootstrap_samples = 1),
               "bootstrap_samples")
  model_set <- fit_propensity_models(instability_formula, lalonde,
                                     algorithms = "logistic_regression")
  no_bootstrap <- match_each_model(model_set, lalonde, replace = FALSE)
  expect_error(bootstrap_rank_instability(no_bootstrap), "no bootstrap block")
})
