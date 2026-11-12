scenarios <- c("S1", "S2", "S3", "S4", "S5")
expected_columns <- function(scenario) {
  c("id", paste0("x", 1:6), if (scenario == "S5") paste0("b", 1:4),
    "treat", "y", "mu0", "tau", "true_logit", "true_propensity")
}
outcome_free_columns <- function(scenario) {
  c("id", paste0("x", 1:6), if (scenario == "S5") paste0("b", 1:4), "treat", "true_logit",
    "true_propensity")
}

test_that("lambda changes only mu0 and y, and only in scenarios with nonlinear terms", {
  for (scenario in scenarios) {
    lambda_zero <- simulate_example_data(2000, scenario, lambda = 0, seed = 7)
    lambda_harm <- simulate_example_data(2000, scenario, lambda = 1.5, seed = 7)
    free <- outcome_free_columns(scenario)
    expect_identical(lambda_zero[free], lambda_harm[free], info = scenario)
    expect_identical(lambda_zero$tau, lambda_harm$tau, info = scenario)
    if (scenario %in% c("S2", "S3", "S5")) {
      expect_false(isTRUE(all.equal(lambda_zero$mu0, lambda_harm$mu0)), info = scenario)
      expect_false(isTRUE(all.equal(lambda_zero$y, lambda_harm$y)), info = scenario)
    } else {
      expect_identical(lambda_zero$mu0, lambda_harm$mu0, info = scenario)
      expect_identical(lambda_zero$y, lambda_harm$y, info = scenario)
    }
  }
})

test_that("treated share is hit at n = 100000 for every scenario", {
  for (scenario in scenarios) {
    data <- simulate_example_data(100000, scenario, seed = 3)
    expect_equal(mean(data$true_propensity), 0.20, tolerance = 1e-6, info = scenario)
    expect_lt(abs(mean(data$treat) - 0.20), 0.005, label = paste(scenario, "realised share"))
  }
  data_40 <- simulate_example_data(100000, "S4", treated_share = 0.40, seed = 3)
  expect_lt(abs(mean(data_40$treat) - 0.40), 0.005)
})

test_that("same seed reproduces, different seeds differ, caller's RNG state is restored", {
  first <- simulate_example_data(300, "S3", lambda = 1.5, seed = 11)
  second <- simulate_example_data(300, "S3", lambda = 1.5, seed = 11)
  other <- simulate_example_data(300, "S3", lambda = 1.5, seed = 12)
  expect_identical(first, second)
  expect_false(isTRUE(all.equal(first$x1, other$x1)))
  expect_false(isTRUE(all.equal(first$treat, other$treat)))
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(simulate_example_data(50, "S1", seed = 1))
  expect_identical(runif(1), before)
})

test_that("columns, types, attributes and no missing values", {
  for (scenario in scenarios) {
    data <- simulate_example_data(500, scenario, lambda = 1.5, seed = 2)
    expect_s3_class(data, "data.frame")
    expect_identical(names(data), expected_columns(scenario))
    expect_equal(nrow(data), 500)
    expect_false(anyNA(data))
    expect_type(data$id, "integer")
    expect_type(data$treat, "integer")
    expect_true(all(data$treat %in% c(0L, 1L)))
    for (column in c(paste0("x", 1:6), "y", "mu0", "tau", "true_logit", "true_propensity")) {
      expect_type(data[[column]], "double")
    }
    if (scenario == "S5") {
      for (column in paste0("b", 1:4)) expect_true(all(data[[column]] %in% c(0L, 1L)))
      expect_lt(mean(data$b4), 0.12)
    }
    expect_true(all(data$true_propensity > 0 & data$true_propensity < 1))
    expect_equal(attr(data, "scenario"), scenario)
    expect_equal(attr(data, "lambda"), 1.5)
    expect_equal(attr(data, "treated_share"), 0.20)
    expect_equal(attr(data, "seed"), 2)
    expect_false(attr(data, "heterogeneity"))
    expect_type(attr(data, "alpha"), "double")
  }
})

test_that("Scenario 4 has poorer overlap than Scenario 1: larger spread of the true logit", {
  s1 <- simulate_example_data(5000, "S1", seed = 5)
  s4 <- simulate_example_data(5000, "S4", seed = 5)
  expect_gt(sd(s4$true_logit), 2 * sd(s1$true_logit))
})

test_that("heterogeneity changes tau and y only", {
  constant <- simulate_example_data(1000, "S2", lambda = 1.5, seed = 4, heterogeneity = FALSE)
  varying <- simulate_example_data(1000, "S2", lambda = 1.5, seed = 4, heterogeneity = TRUE)
  unchanged <- setdiff(names(constant), c("tau", "y"))
  expect_identical(constant[unchanged], varying[unchanged])
  expect_true(all(constant$tau == 1))
  expect_equal(varying$tau, 1 + 0.8 * varying$x1)
  expect_false(isTRUE(all.equal(constant$y, varying$y)))
  expect_equal(varying$y - constant$y, (varying$tau - 1) * varying$treat)
})

test_that("bad inputs give clear errors", {
  expect_error(simulate_example_data(1, "S1"), "`n`")
  expect_error(simulate_example_data(100, "S9"), "should be one of")
  expect_error(simulate_example_data(100, "S1", lambda = "a"), "`lambda`")
  expect_error(simulate_example_data(100, "S1", treated_share = 1.2), "`treated_share`")
  expect_error(simulate_example_data(100, "S1", seed = 1.5), "`seed`")
  expect_error(simulate_example_data(100, "S1", heterogeneity = NA), "`heterogeneity`")
})
