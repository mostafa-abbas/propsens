# Bootstrap rank instability of the logistic propensity model

Descriptive. The sampling noise of the analyst's logistic model:
logistic regression with the formula as written (including any
interactions or transformations) is refitted on `bootstrap_samples`
nonparametric bootstrap samples of all patients, the original patients
are ranked by each refit's logit (percentile rank
`rank(logit, ties = "average") / n`, on a 0 to 1 scale), and for each
patient the standard deviation of the rank over the refits is taken. The
study-level value is the mean over treated patients. A bootstrap sample
with fewer than 10 treated or 10 untreated patients is skipped (the rule
of
[`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md)).

## Usage

``` r
bootstrap_rank_instability(
  x,
  data = NULL,
  bootstrap_samples = 20,
  cores = 1,
  seed = 1,
  clip_bootstrap = FALSE
)

# S3 method for class 'bootstrap_rank_instability'
print(x, ...)
```

## Arguments

- x:

  One of: a `propensity_sensitivity` object or a `match_set` from
  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md)
  (the stored bootstrap ranks are used, nothing is refitted); a
  `matchit` object (its formula is used); or a two-sided formula such as
  `treat ~ x1 + x2`.

- data:

  The data frame. Required with a formula; with a `matchit` object it is
  looked up from the object's call when `NULL`. Not used for stored
  ranks.

- bootstrap_samples:

  Number of bootstrap refits (default 20). Not used for stored ranks.

- cores:

  Number of cores (default 1). Not used for stored ranks.

- seed:

  Integer seed from which every other seed is derived; the same seed
  gives the bootstrap samples of
  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md).
  Not used for stored ranks.

- clip_bootstrap:

  `FALSE` (default) or `TRUE`: clip the refits' scores to
  `[0.001, 0.999]`, the default `clip` of
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md).
  Not used for stored ranks.

- ...:

  Not used.

## Value

An object of class `"bootstrap_rank_instability"`: a list with
`patients` (one row per patient: `id`, `treated`, `rank_instability`),
`summary` (`mean_treated`, `samples_used`, `samples_requested`,
`samples_skipped`) and `settings`. `mean_treated` is `NA` when fewer
than two bootstrap samples are usable.

## Details

It uses the covariates and the treatment only, no outcome. Report it
next to rank disagreement: it shows how much the logistic model's ranks
move from sampling alone. It describes the analyst's logistic model
only, whatever distance the `matchit` object used. The refits' scores
are not clipped by default, unlike the cross-fitted scores (see
`clip_bootstrap`). These diagnostics have been studied in the author's
simulations and benchmark analyses; a methods paper is in preparation.
No peer-reviewed evaluation exists yet.

## Examples

``` r
data("lalonde", package = "MatchIt")
instability <- bootstrap_rank_instability(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, bootstrap_samples = 5
)
instability
#> Bootstrap rank instability of the logistic model (mean over treated): 0.03605
#> SD of the percentile rank over 5 bootstrap refits
instability$summary$mean_treated
#> [1] 0.0360503
```
