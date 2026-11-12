# Per-patient rank disagreement across propensity models

For each patient, `d` is the standard deviation over the kept models of
the percentile rank `r_m(i) = rank(logit_m(i)) / n` (0 to 1 scale). Rank
disagreement, the study mean of `d` over treated patients, is the
primary summary; it is a study-level summary, not a per-patient flag.

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
`logit_disagreement`), `summary` (`mean_d_treated`,
`mean_d_between_treated`, `mean_d_all`,
`median_between_within_ratio_treated`, `mean_d_logit_treated`),
`models_used`, `quality_gate_left_fewer_than_two_models` and `settings`.

## Details

Two exploratory companions, which are not a decomposition of `d` (`d^2`
is not `d_between^2 + d_within^2`): `d_between` is the standard
deviation over algorithms of each algorithm's mean rank, so that every
algorithm counts once however many seeds it has; its study mean over
treated patients is reported as "rank disagreement, equal weight per
algorithm" (`NA` with one algorithm). `d_within` is the square root of
the mean, over algorithms with at least two seeds, of the variance over
seeds (`NA` when no algorithm has two seeds). The ratio
`d_between / d_within` is the disagreement between algorithms divided by
the disagreement between seeds of the same algorithm. With three seeds
per algorithm, seed noise alone gives a ratio of about 0.4; larger
values point to disagreement between algorithms beyond seed noise.

`d` is a standard deviation over models, so it depends on how many seeds
each algorithm contributes: with the same algorithms, more seeds of one
algorithm pull `d` towards that algorithm's seed noise. Values are
comparable only between runs with the same model set and the same kept
algorithms. `d_logit`, the standard deviation of the logits, is a
logit-scale version; calibration differences inflate it without changing
who is near whom. Standard deviations use the n - 1 denominator.

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
  data = lalonde, seeds = 1:2, n_trees = 100,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
)
disagreement <- rank_disagreement(model_set)
disagreement
#> Rank disagreement: 4 model(s) in use (after the log-loss gate)
#> Mean d over treated (primary): 0.02894; over all patients: 0.0445
#> Equal weight per algorithm, mean d_between over treated (exploratory): 0.03256
#> Median d_between / d_within over treated (exploratory): 3.522
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
