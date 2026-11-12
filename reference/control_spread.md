# Per-patient control spread across propensity models

Descriptive. ("Control" in the name means the untreated partner.) For
each treated patient, collect the distinct untreated partners matched to
it by the kept models and compute `s`, the mean pairwise Mahalanobis
distance between their covariate vectors, with the covariance matrix of
all n patients. `s = 0` when every model picks the same untreated
partner, and `NA` when the patient is unmatched in every model.
`s_scaled = s / D`, where `D` is the mean over matched treated patients
of the Mahalanobis distance between a treated patient and its
logistic-regression match (the first model in use when the match set has
no logistic regression model): 1 means the candidate partners differ
from each other as much as a typical match differs from its treated
patient. The scale model is the logistic regression match even when
logistic regression is not one of the kept models;
[`print()`](https://rdrr.io/r/base/print.html) names it. When `D` is 0
(every scale pair is identical, for example with a very tight caliper)
`s_scaled` is `NA`, not infinite. The oracle entry and the bootstrap
block are never used.

## Usage

``` r
control_spread(match_set, covariates, use_gate = TRUE, models = NULL)

# S3 method for class 'control_spread'
print(x, ...)
```

## Arguments

- match_set:

  A `match_set` from
  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md).

- covariates:

  A data frame with one row per patient, in the data's row order,
  holding the covariates that define the distance (usually those of the
  propensity formula).

- use_gate:

  `TRUE` (default) to use the models that passed the log-loss gate;
  `FALSE` to use every model.

- models:

  Optional character vector of model ids to use instead of a gate; it
  overrides `use_gate`. The wrapper
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  passes the final quality-gate set (log-loss and balance) this way.

- x:

  A `control_spread` object.

- ...:

  Not used.

## Value

An object of class `"control_spread"`: a list with `patients` (one row
per patient: `row`, `id`, `treated`, `distinct_controls`,
`control_spread`, `control_spread_scaled`; `NA` for untreated patients),
`summary` (`mean_scaled_treated`, `p90_scaled_treated`,
`scale_distance`, `scale_model`), `models_used`,
`quality_gate_left_fewer_than_two_models` and `settings` (with the
dropped columns).

## Covariates

Factors and logical columns are expanded to dummies (treatment
contrasts). Character columns are an error: convert them to factors.
Constant columns and linearly dependent columns (found with a pivoted QR
decomposition of the standardized design) are dropped with a message.
The Mahalanobis distance is `sqrt((a - b)' S^-1 (a - b))`, the square
root of what
[`stats::mahalanobis()`](https://rdrr.io/r/stats/mahalanobis.html)
returns.

## Quality gate

"Kept models" in this help page are the models the function uses: the
models given in `models` (as
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
passes the models kept by its quality gate), or with `use_gate = TRUE`
the models that passed the log-loss gate, or every model. If fewer than
two passed, all models are used and
`quality_gate_left_fewer_than_two_models` is `TRUE`. The flag reports
the gate's outcome whatever `use_gate` is.

## Examples

``` r
data("lalonde", package = "MatchIt")
model_set <- fit_propensity_models(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, seeds = 1, n_trees = 100,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
)
match_set <- match_each_model(model_set, lalonde, replace = TRUE)
covariates <- lalonde[, c("age", "educ", "race", "married", "nodegree", "re74", "re75")]
spread <- control_spread(match_set, covariates)
spread
#> Control spread: 2 model(s) in use (after the log-loss gate)
#> Scaled control spread over treated (descriptive): mean 1.035, 90th percentile 1.927
#> Scale: mean Mahalanobis distance treated to match under logistic_regression_1 = 2.584
head(spread$patients[spread$patients$treated == 1, ])
#>   row   id treated distinct_controls control_spread control_spread_scaled
#> 1   1 NSW1       1                 2      4.3937831             1.7004394
#> 2   2 NSW2       1                 2      5.8022919             2.2455468
#> 3   3 NSW3       1                 1      0.0000000             0.0000000
#> 4   4 NSW4       1                 2      2.4397568             0.9442110
#> 5   5 NSW5       1                 2      0.9652184             0.3735495
#> 6   6 NSW6       1                 2      2.1516074             0.8326943
```
