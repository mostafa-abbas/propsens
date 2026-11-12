data("lalonde", package = "MatchIt")
plot_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75
plot_result <- propensity_sensitivity(plot_formula, data = lalonde,
                                      algorithms = c("logistic_regression", "random_forest",
                                                     "gradient_boosting"),
                                      replace = TRUE, seeds = 1:2, n_trees = 50,
                                      bootstrap_samples = 5)

layer_geoms <- function(plot) {
  vapply(plot$layers, function(layer) class(layer$geom)[1], character(1))
}

test_that("each type returns a ggplot object with the expected geoms", {
  patients <- plot(plot_result, type = "patients")
  models <- plot(plot_result, type = "models")
  recurrence <- plot(plot_result, type = "recurrence")
  for (p in list(patients, models, recurrence)) expect_s3_class(p, "ggplot")
  expect_setequal(layer_geoms(patients), c("GeomVline", "GeomHline", "GeomPoint"))
  expect_setequal(layer_geoms(models), c("GeomVline", "GeomPoint"))
  expect_setequal(layer_geoms(recurrence), c("GeomBar", "GeomVline"))
  expect_error(plot(plot_result, type = "pie"), "should be one of")
})

test_that("the patients plot has one point per matched treated patient", {
  built <- ggplot2::ggplot_build(plot(plot_result, type = "patients"))
  points <- built$data[[3]]
  treated <- plot_result$patients[plot_result$patients$treated == 1, ]
  expect_equal(nrow(points), sum(!is.na(treated$control_spread_scaled)))
  expect_equal(built$data[[1]]$xintercept[1], stats::median(treated$rank_disagreement))
})

test_that("the models plot marks the reference and gate failures", {
  built <- ggplot2::ggplot_build(plot(plot_result, type = "models"))
  points <- built$data[[2]]
  expect_equal(nrow(points), 4 * nrow(plot_result$models))
  passed <- plot_result$models$passed_quality_gate %in% TRUE
  # shape 1 = hollow (failed), 16 = solid (passed)
  expect_setequal(unique(points$shape), unique(ifelse(passed, 16, 1)))
  expect_equal(built$data[[1]]$xintercept, plot_result$settings$balance_limit)
})

test_that("no error with NA rows (unmatched treated, missing support)", {
  with_na <- plot_result
  treated_rows <- which(with_na$patients$treated == 1)
  with_na$patients$control_spread_scaled[treated_rows[1:5]] <- NA
  with_na$patients$outside_support[treated_rows] <- NA
  with_na$patients$pair_recurrence_bootstrap[treated_rows[1:5]] <- NA
  with_na$models$smd_max_extended[1] <- NA
  for (type in c("patients", "models", "recurrence")) {
    expect_no_error(ggplot2::ggplot_build(plot(with_na, type = type)))
  }
  built <- ggplot2::ggplot_build(plot(with_na, type = "patients"))
  matched_treated <- !is.na(with_na$patients$control_spread_scaled[treated_rows])
  expect_equal(nrow(built$data[[3]]), sum(matched_treated))
})

test_that("snapshot of the built plot data (not the image)", {
  round_columns <- function(frame, columns) {
    frame <- frame[, intersect(columns, names(frame)), drop = FALSE]
    for (column in names(frame)) {
      if (is.numeric(frame[[column]])) frame[[column]] <- signif(frame[[column]], 4)
    }
    frame
  }
  built_patients <- ggplot2::ggplot_build(plot(plot_result, type = "patients"))
  built_models <- ggplot2::ggplot_build(plot(plot_result, type = "models"))
  built_recurrence <- ggplot2::ggplot_build(plot(plot_result, type = "recurrence"))
  expect_snapshot_value(list(
    patients = round_columns(head(built_patients$data[[3]], 5), c("x", "y", "colour")),
    models = round_columns(built_models$data[[2]], c("x", "y", "colour", "shape", "PANEL")),
    recurrence = round_columns(built_recurrence$data[[1]], c("x", "count", "PANEL")),
    recurrence_lines = round_columns(built_recurrence$data[[2]], c("xintercept", "PANEL"))
  ), style = "deparse")
})
