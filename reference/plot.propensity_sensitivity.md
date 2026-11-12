# Plots of a propensity model sensitivity result

Three plots, built with ggplot2 and a colour-blind-safe (Okabe-Ito)
palette, readable at 8 cm width.

## Usage

``` r
# S3 method for class 'propensity_sensitivity'
plot(x, type = c("patients", "models", "recurrence"), ...)
```

## Arguments

- x:

  A `propensity_sensitivity` object.

- type:

  One of `"patients"`, `"models"`, `"recurrence"`.

- ...:

  Not used.

## Value

A ggplot object.

## Details

- `type = "patients"`: how the disagreement is spread over the treated
  patients: one point per treated patient; x = rank disagreement, y =
  scaled control spread, colour = outside the support of the untreated
  patients' scores; dashed lines at the study medians. It is a
  description of the spread, not a set of flags for individual patients.
  Treated patients unmatched in every kept model (caliper) are left out.

- `type = "models"`: the four balance statistics by model, one panel per
  statistic, each axis starting at 0. Kept models are filled circles,
  the others hollow; your own model (from a `matchit` object) is a
  diamond in its own colour, filled when its balance is within the limit
  (it is judged on balance only and is never kept); the dashed line in
  the main-terms panel is the balance limit of the quality gate.

- `type = "recurrence"`: histograms of pair recurrence over treated
  patients, across the kept models (all models after a fallback; floor 1
  / number of models and ceiling 1 drawn; pair instability is 1 minus
  its mean, and it is high whenever many untreated patients have similar
  scores, even when the models agree) and across bootstrap refits of the
  logistic model with the corrected denominator (ceiling 1).

## Examples

``` r
data("lalonde", package = "MatchIt")
matched <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                            data = lalonde, distance = "glm", link = "linear.logit",
                            replace = TRUE)
sensitivity <- propensity_sensitivity(
  matched, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
  seeds = 1, n_trees = 100, bootstrap_samples = 5
)
#> Warning: Fewer than two models passed the quality gate, so there is no set of similarly good models to compare. Every model is used below; read the values as a description of all models, not of a kept set.
plot(sensitivity, type = "patients")

# \donttest{
# Each further plot adds a few seconds of drawing time.
plot(sensitivity, type = "models")

plot(sensitivity, type = "recurrence")

# }
```
