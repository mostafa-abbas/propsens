# Repeat the diagnostic with several values of `seed`

The fold split and every tree model of
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
depend on `seed`, and on small data the kept models and the summaries
can change noticeably from one seed to the next. `seed_sensitivity()`
reruns
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
once per value of `seeds` (here: values of its `seed` argument) with
every other argument unchanged, and returns one row per seed with the
models kept, whether the gate fell back, and the main summaries, plus
their range. Report the whole table.

## Usage

``` r
seed_sensitivity(x, seeds = 1:3, ..., model_seeds = NULL, keep_results = FALSE)

# S3 method for class 'seed_sensitivity'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `matchit` object or a two-sided formula, as for
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md).

- seeds:

  Values of `seed` to run (default `1:3`). Each value is one full run.

- ...:

  Further arguments for
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  (for example `data`, `algorithms`, `replace`, `caliper`, `cores`),
  except `seed`.

- model_seeds:

  The `seeds` argument of
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  (the number of seeds per tree algorithm, as `1:k`); `NULL` (default)
  keeps its default `1:3`.

- keep_results:

  `TRUE` to keep every full result in `results` (large); default
  `FALSE`.

- digits:

  Significant digits for printed numbers.

## Value

An object of class `"seed_sensitivity"`: a list with `table` (one row
per seed: `seed`, `models_kept`, `models_fitted`, `algorithms_kept`,
`kept_by_algorithm`, `fallback`, `rank_disagreement`,
`rank_disagreement_equal_weight`, `pair_instability`,
`bootstrap_rank_instability`), `range` (minimum and maximum of each
summary over the seeds) and `results` (`NULL` unless
`keep_results = TRUE`).

## Details

The R warnings of each run about the gate fallback and about a single
kept algorithm are not repeated; they are columns of the table
(`fallback`, `algorithms_kept`).

## Examples

``` r
data("lalonde", package = "MatchIt")
by_seed <- seed_sensitivity(
  treat ~ age + educ + race + married + nodegree + re74 + re75, seeds = 1:2,
  data = lalonde, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
  model_seeds = 1, n_trees = 50, bootstrap_samples = 2
)
by_seed
#> Propensity model sensitivity over 2 values of `seed`
#>  seed models_kept models_fitted algorithms_kept
#>     1           3             3               3
#>     2           3             3               3
#>                                            kept_by_algorithm fallback
#>  logistic regression 1, random forest 1, gradient boosting 1     TRUE
#>  logistic regression 1, random forest 1, gradient boosting 1     TRUE
#>  rank_disagreement rank_disagreement_equal_weight pair_instability
#>             0.0692                         0.0692            0.659
#>             0.0677                         0.0677            0.659
#>  bootstrap_rank_instability
#>                      0.0230
#>                      0.0348
#> 
#> Range over the seeds:
#>                         summary minimum maximum
#>               rank_disagreement  0.0677  0.0692
#>  rank_disagreement_equal_weight  0.0677  0.0692
#>                pair_instability  0.6590  0.6590
#>      bootstrap_rank_instability  0.0230  0.0348
#> Some seeds fell back to every model or kept a single algorithm; their values
#>    do not describe a set of similarly good models from several algorithms.
#> Report every row; the values depend on the seed as much as shown here.
```
