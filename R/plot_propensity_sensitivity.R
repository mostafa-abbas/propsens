#' Plots of a propensity model sensitivity result
#'
#' Three plots, built with ggplot2 and a colour-blind-safe
#' (Okabe-Ito) palette, readable at 8 cm width.
#'
#' - `type = "patients"`: one point per treated patient; x = rank disagreement, y = scaled
#'   control spread, colour = outside the support of the untreated patients' scores; dashed
#'   lines at the study medians. Upper right = the models disagree about this patient and pick
#'   different kinds of untreated partners. Treated patients unmatched in every model (caliper)
#'   are left out.
#' - `type = "models"`: the four balance statistics by model, one panel per
#'   statistic. Models that fail the quality gate are hollow; the analyst's reference model is
#'   drawn in its own colour; the dashed line in the main-terms panel is the balance limit of
#'   the quality gate.
#' - `type = "recurrence"`: histograms of pair recurrence over treated patients, across models
#'   (floor 1 / number of models and ceiling 1 drawn) and across bootstrap resamples with the
#'   corrected denominator (ceiling 1; the expected ceiling 1 - e^(-1) about 0.63 of the uncorrected
#'   version is drawn for reference).
#'
#' @param x A `propensity_sensitivity` object.
#' @param type One of `"patients"`, `"models"`, `"recurrence"`.
#' @param ... Not used.
#' @return A ggplot object.
#' @examples
#' data("lalonde", package = "MatchIt")
#' matched <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
#'                             data = lalonde, distance = "glm", link = "linear.logit",
#'                             replace = TRUE)
#' sensitivity <- propensity_sensitivity(
#'   matched, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
#'   seeds = 1, n_trees = 100, bootstrap_samples = 5
#' )
#' plot(sensitivity, type = "patients")
#' \donttest{
#' # Each further plot adds a few seconds of drawing time.
#' plot(sensitivity, type = "models")
#' plot(sensitivity, type = "recurrence")
#' }
#' @importFrom ggplot2 .data
#' @export
plot.propensity_sensitivity <- function(x, type = c("patients", "models", "recurrence"), ...) {
  type <- match.arg(type)
  switch(type,
         patients = plot_patients(x),
         models = plot_models(x),
         recurrence = plot_recurrence(x))
}

# Okabe-Ito colours: blue, vermillion, grey, bluish green.
okabe_ito <- c(blue = "#0072B2", vermillion = "#D55E00", grey = "#999999",
               green = "#009E73")

plot_theme <- function() {
  ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 7),
                   panel.grid.minor = ggplot2::element_blank())
}

plot_patients <- function(x) {
  patients <- x$patients[x$patients$treated == 1, , drop = FALSE]
  patients <- patients[!is.na(patients$rank_disagreement) &
                         !is.na(patients$control_spread_scaled), , drop = FALSE]
  patients$support <- factor(ifelse(is.na(patients$outside_support), "unknown",
                                    ifelse(patients$outside_support, "yes", "no")),
                             levels = c("no", "yes", "unknown"))
  medians <- data.frame(
    rank_disagreement = stats::median(patients$rank_disagreement),
    control_spread_scaled = stats::median(patients$control_spread_scaled)
  )
  ggplot2::ggplot(patients, ggplot2::aes(x = .data$rank_disagreement,
                                         y = .data$control_spread_scaled)) +
    ggplot2::geom_vline(data = medians, ggplot2::aes(xintercept = .data$rank_disagreement),
                        linetype = "dashed", colour = okabe_ito[["grey"]]) +
    ggplot2::geom_hline(data = medians, ggplot2::aes(yintercept = .data$control_spread_scaled),
                        linetype = "dashed", colour = okabe_ito[["grey"]]) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$support), size = 1.2, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c(no = okabe_ito[["blue"]],
                                            yes = okabe_ito[["vermillion"]],
                                            unknown = okabe_ito[["grey"]]),
                                 drop = TRUE, name = "Outside support") +
    ggplot2::labs(x = "Rank disagreement across models",
                  y = "Scaled control spread",
                  caption = paste0("One point per treated patient; dashed lines: medians.\n",
                                   "Upper right = the models disagree about this patient\n",
                                   "and pick different kinds of controls.")) +
    plot_theme()
}

plot_models <- function(x) {
  statistics <- c(smd_max_main = "Max |SMD| main",
                  smd_max_extended = "Max |SMD| extended",
                  energy_distance = "Energy distance",
                  ks_max = "Max KS")
  models <- x$models
  long <- do.call(rbind, lapply(names(statistics), function(statistic) {
    data.frame(model_id = models$model_id, statistic = statistics[[statistic]],
               value = models[[statistic]],
               kind = ifelse(models$is_reference, "reference (analyst)", "refitted model"),
               gate = ifelse(models$passed_quality_gate %in% TRUE, "passed", "failed"),
               stringsAsFactors = FALSE)
  }))
  long$statistic <- factor(long$statistic, levels = unname(statistics))
  long$model_id <- factor(long$model_id, levels = rev(models$model_id))
  long$kind <- factor(long$kind, levels = c("refitted model", "reference (analyst)"))
  long$gate <- factor(long$gate, levels = c("passed", "failed"))
  threshold <- data.frame(statistic = factor(statistics[["smd_max_main"]],
                                             levels = unname(statistics)),
                          value = x$settings$balance_limit)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$model_id)) +
    ggplot2::geom_vline(data = threshold, ggplot2::aes(xintercept = .data$value),
                        linetype = "dashed", colour = okabe_ito[["grey"]]) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$kind, shape = .data$gate), size = 2,
                        na.rm = TRUE) +
    ggplot2::scale_shape_manual(values = c(passed = 16, failed = 1), drop = FALSE,
                                name = "Quality gate") +
    ggplot2::scale_colour_manual(values = c("refitted model" = okabe_ito[["blue"]],
                                            "reference (analyst)" = okabe_ito[["vermillion"]]),
                                 drop = TRUE, name = NULL) +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::facet_wrap(~ statistic, scales = "free_x", ncol = 2) +
    ggplot2::labs(x = "Value after matching (lower = better)", y = NULL,
                  caption = paste0("Extended: with squares, interactions.\n",
                                   "Hollow = fails the quality gate;\n",
                                   "dashed line = balance limit of the quality gate.")) +
    plot_theme()
}

plot_recurrence <- function(x) {
  treated <- x$patients[x$patients$treated == 1, , drop = FALSE]
  panels <- c(models = "Across models", bootstrap = "Across bootstrap resamples")
  long <- rbind(
    data.frame(panel = panels[["models"]], recurrence = treated$pair_recurrence_models),
    data.frame(panel = panels[["bootstrap"]], recurrence = treated$pair_recurrence_bootstrap)
  )
  long <- long[!is.na(long$recurrence), , drop = FALSE]
  long$panel <- factor(long$panel, levels = unname(panels))
  lines <- data.frame(
    panel = factor(c(panels[["models"]], panels[["models"]], panels[["bootstrap"]],
                     panels[["bootstrap"]]), levels = unname(panels)),
    value = c(x$summary$pair_recurrence_floor_models, 1,
              x$summary$pair_recurrence_expected_ceiling_bootstrap_uncorrected, 1),
    line = c("floor (1 / models)", "ceiling", "uncorrected ceiling (0.63)", "ceiling")
  )
  ggplot2::ggplot(long, ggplot2::aes(x = .data$recurrence)) +
    ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, fill = okabe_ito[["blue"]],
                            colour = "white") +
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = .data$value,
                                                   linetype = .data$line),
                        colour = okabe_ito[["vermillion"]]) +
    ggplot2::scale_linetype_manual(values = c("floor (1 / models)" = "dotted",
                                              "ceiling" = "solid",
                                              "uncorrected ceiling (0.63)" = "dashed"),
                                   name = NULL) +
    ggplot2::guides(linetype = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y", drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, 1.05)) +
    ggplot2::labs(x = "Pair recurrence", y = "Treated patients",
                  caption = paste0("Share of runs that choose the modal control. ",
                                   "Descriptive index:\n",
                                   "it measures crowding of controls, not model sensitivity.")) +
    plot_theme()
}
