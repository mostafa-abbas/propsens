# Per-patient rank disagreement across propensity models

For each patient, `d` is the standard deviation over the models in use
of the percentile rank `r_m(i) = rank(logit_m(i)) / n`. The
decomposition separates disagreement between algorithms from seed noise
within an algorithm: `d_between` is the standard deviation over
algorithms of the per-algorithm mean rank, and `d_within` is the square
root of the mean over algorithms with at least two seeds of the variance
over seeds (`NA` when no algorithm has two seeds). A ratio
`d_between / d_within` near 1 means the algorithms disagree with
themselves as much as with each other: noise, not model dependence.
`d_logit`, the standard deviation of the logits, is the secondary
logit-scale version; calibration differences inflate it without changing
who is near whom. Standard deviations use the n - 1 denominator.

## Usage

``` r
rank_disagreement(model_set, use_gate = TRUE, models = NULL)

# S3 method for class 'rank_disagreement'
print(x, ...)
```

## Arguments

- model_set:

  A `propensity_model_set` from
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md).

- use_gate:

  `TRUE` (default) to use the models that passed the log-loss gate;
  `FALSE` to use every model.

- models:

  Optional character vector of model ids to use instead of a gate; it
  overrides `use_gate`. The wrapper
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  passes the final quality-gate set (log-loss and balance) this way.

- x:

  A `rank_disagreement` object.

- ...:

  Not used.

## Value

An object of class `"rank_disagreement"`: a list with `patients` (one
row per patient: `row`, `treated`, `rank_disagreement`,
`rank_disagreement_between`, `rank_disagreement_within`,
`logit_disagreement`), `summary` (`mean_d_treated`, `mean_d_all`,
`median_between_within_ratio_treated`, `mean_d_logit_treated`),
`models_used`, `quality_gate_left_fewer_than_two_models` and `settings`.

## Quality gate

With `use_gate = TRUE` only the models that passed the log-loss gate are
used. If fewer than two passed, all models are used and
`quality_gate_left_fewer_than_two_models` is `TRUE` . The flag reports
the gate's outcome whatever `use_gate` is.

## Examples

``` r
data("lalonde", package = "MatchIt")
model_set <- fit_propensity_models(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, seeds = 1:2, n_trees = 100,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
)
disagreement <- rank_disagreement(model_set)
disagreement
#> Rank disagreement: 4 model(s) in use (after the log-loss gate)
#> Mean d over treated: 0.02894; over all patients: 0.0445
#> Median d_between / d_within over treated: 3.522 (near 1 = seed noise, not model dependence)
#> Mean d_logit over treated: 0.4483
head(disagreement$patients)
#>   row treated rank_disagreement rank_disagreement_between
#> 1   1       1       0.008410387               0.000000000
#> 2   2       1       0.033546140               0.040883210
#> 3   3       1       0.020904719               0.024760254
#> 4   4       1       0.002338991               0.002015369
#> 5   5       1       0.018370669               0.021593245
#> 6   6       1       0.013287597               0.015835046
#>   rank_disagreement_within logit_disagreement
#> 1              0.010300579          0.1534271
#> 2              0.004071661          0.1661779
#> 3              0.006514658          0.7622771
#> 4              0.002035831          1.4968027
#> 5              0.006320918          0.2557289
#> 6              0.003753886          0.4665527
```
