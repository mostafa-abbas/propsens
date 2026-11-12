# Covariates of the balance check, the quality gate and the control spread: the right-hand
# side is resolved with R's formula rules, so "." respects exclusions.

make_study <- function(n = 300, seed = 11) {
  withr::with_seed(seed, {
    age <- round(stats::rnorm(n, 60, 10))
    diabetes <- stats::rbinom(n, 1, 0.3)
    sex <- factor(sample(c("female", "male"), n, replace = TRUE))
    treat <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.04 * (age - 60) + 0.6 * diabetes))
    # An outcome strongly tied to treatment: it would dominate the balance check if used.
    death <- stats::rbinom(n, 1, ifelse(treat == 1, 0.8, 0.1))
    data.frame(treat = treat, death = death, patient_id = seq_len(n), age = age,
               diabetes = diabetes, sex = sex)
  })
}

covariate_names <- function(formula, data) {
  names(select_covariates(NULL, formula, data))
}

test_that("'.' with exclusions keeps only the remaining columns (formula and matchit object)", {
  study <- make_study()
  formula <- treat ~ . - death - patient_id
  expect_identical(covariate_names(formula, study), c("age", "diabetes", "sex"))
  matched <- MatchIt::matchit(formula, data = study)
  design <- read_matchit_design(matched, study, environment())
  expect_identical(covariate_names(design$formula, study), c("age", "diabetes", "sex"))
})

test_that("'.' alone gives every column except the treatment", {
  study <- make_study()
  expect_identical(covariate_names(treat ~ ., study), setdiff(names(study), "treat"))
})

test_that("transformations and interactions give their source columns", {
  study <- make_study()
  expect_identical(covariate_names(treat ~ age + I(age^2), study), "age")
  expect_identical(covariate_names(treat ~ age + diabetes + age:diabetes, study),
                   c("age", "diabetes"))
  expect_identical(covariate_names(treat ~ factor(diabetes) + sex, study), c("diabetes", "sex"))
  expect_identical(covariate_names(treat ~ splines::ns(age, df = 3) + sex, study),
                   c("age", "sex"))
})

test_that("an excluded outcome never enters the balance check (end to end)", {
  study <- make_study()
  run <- function(x, data = study) {
    propensity_sensitivity(x, data = data, algorithms = c("logistic_regression",
                                                          "gradient_boosting"),
                           seeds = 1, n_trees = 50, bootstrap_samples = 0, seed = 2)
  }
  dot <- run(treat ~ . - death - patient_id)
  explicit <- run(treat ~ age + diabetes + sex)
  expect_identical(dot$settings$covariate_names, c("age", "diabetes", "sex"))
  expect_false("death" %in% dot$settings$covariate_names)
  expect_equal(dot$models$smd_max_main, explicit$models$smd_max_main)
  expect_identical(dot$models$passed_quality_gate, explicit$models$passed_quality_gate)
  expect_equal(dot$summary$table, explicit$summary$table)

  matched <- MatchIt::matchit(treat ~ . - death - patient_id, data = study)
  from_matchit <- run(matched)
  expect_identical(from_matchit$settings$covariate_names, c("age", "diabetes", "sex"))
})

test_that("custom covariates with differing row names stop with a clear error", {
  study <- make_study()
  rownames(study) <- paste0("patient_", seq_len(nrow(study)))
  covariates <- study[, c("age", "diabetes", "sex")]
  expect_identical(select_covariates(covariates, treat ~ age, study), covariates)
  shuffled <- withr::with_seed(3, covariates[sample(nrow(covariates)), ])
  expect_error(select_covariates(shuffled, treat ~ age, study),
               "row names of `covariates` differ from those of `data`")
  expect_error(propensity_sensitivity(treat ~ age + diabetes + sex, data = study,
                                      algorithms = "logistic_regression",
                                      covariates = shuffled),
               "same patients in the same order")
})
