# Study-level balance spread across propensity models

Descriptive. For each model, the model's MatchIt weights from
[`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md)
are passed to cobalt and four balance statistics are computed:

- `smd_max_main`: maximum absolute standardized mean difference over the
  main terms;

- `smd_max_extended`: the same over main terms plus squares and 2-way
  interactions
  (`cobalt::bal.tab(int = interactions, poly = polynomial)`);

- `energy_distance`:
  [`cobalt::bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.html)
  with `stat = "energy.dist"`;

- `ks_max`:
  [`cobalt::bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.html)
  with `stat = "ks.max"`.

## Usage

``` r
balance_spread(
  match_set,
  treatment,
  covariates,
  interactions = TRUE,
  polynomial = 2,
  use_gate = TRUE,
  models = NULL
)

# S3 method for class 'balance_spread'
print(x, ...)
```

## Arguments

- match_set:

  A `match_set` from
  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md).

- treatment:

  The 0/1 treatment vector, in the data's row order; must equal the
  treatment stored in `match_set`.

- covariates:

  A data frame with one row per patient, in the data's row order.
  Factors are split into dummies by cobalt; character columns are an
  error.

- interactions:

  `TRUE` (default) to include all 2-way interactions in
  `smd_max_extended`.

- polynomial:

  Highest power of the continuous covariates in `smd_max_extended`
  (default 2, squares).

- use_gate:

  `TRUE` (default) to compute the spread over the models that passed the
  log-loss gate; `FALSE` to use every model.

- models:

  Optional character vector of model ids to use instead of a gate; it
  overrides `use_gate`. The wrapper
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  passes the final quality-gate set (log-loss and balance) this way.

- x:

  A `balance_spread` object.

- ...:

  Not used.

## Value

An object of class `"balance_spread"`: a list with `models` (one row per
model: `model_id`, `in_use`, `smd_max_main`, `smd_max_extended`,
`energy_distance`, `ks_max`), `spread` (one row per statistic:
`statistic`, `range`, `sd`, `worst_value`, `worst_model`), `summary`,
`models_used`, `quality_gate_left_fewer_than_two_models` and `settings`.

## Details

With many covariates `smd_max_extended` is usually driven by sparse
squares and interactions (a rare level times another covariate) and can
be large even when every main term is well balanced; report the
main-terms range (`smd_max_main`) as the headline balance spread. With a
caliper the standardized mean differences compare the matched treated
patients with their matched untreated partners, while the energy
distance (`estimand = "ATT"`) also measures the distance to the whole
treated group, so it also grows when treated patients are dropped.

The spread of each statistic over the kept models is its range (max -
min), its standard deviation, and the worst (largest) value with the
model that produced it. The table of models lists every model of the
match set with an `in_use` column, so that the balance part of the
quality gate can read every model; the spread uses the kept models only.
The oracle entry and the bootstrap block are never used.

## cobalt settings (verified against cobalt 5.0.0)

`estimand = "ATT"`, `method = "matching"`, `s.d.denom = "treated"`: the
standardization uses the standard deviation of the treated group in the
unweighted sample. `binary = "std"` and `continuous = "std"`: binary
covariates are standardized as well (cobalt's `bal.tab()` default for
binary covariates is the raw difference in proportions, while
`bal.compute(stat = "smd.max")` standardizes them; the package uses the
standardized version throughout so the two agree). Example: on lalonde
after the logistic match with replacement, `married` has a raw
difference in proportions of 0.059 and a standardized difference of
0.152, so the choice decides whether a 0.10 balance threshold is met.
Weighted means, the energy distance and the KS statistics do not depend
on how the weights of the untreated patients are scaled.

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
balance <- balance_spread(match_set, lalonde$treat, covariates)
balance
#> Balance spread: 2 model(s) in use (after the log-loss gate)
#>               model_id in_use smd_max_main smd_max_extended energy_distance
#>  logistic_regression_1  FALSE       0.2710           0.4588         0.09529
#>        random_forest_1   TRUE       0.1600           0.2442         0.05546
#>    gradient_boosting_1   TRUE       0.2854           0.3769         0.08851
#>  ks_max
#>  0.3568
#>  0.1405
#>  0.1351
#> Spread over the models in use:
#>         statistic    range       sd worst_value         worst_model
#>      smd_max_main 0.125300 0.088630     0.28540 gradient_boosting_1
#>  smd_max_extended 0.132700 0.093840     0.37690 gradient_boosting_1
#>   energy_distance 0.033050 0.023370     0.08851 gradient_boosting_1
#>            ks_max 0.005405 0.003822     0.14050     random_forest_1
```
