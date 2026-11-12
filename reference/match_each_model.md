# Match the treated patients once per propensity model

Runs 1:1 nearest-neighbour matching for the average treatment effect in
the treated (ATT) once per model of a
[`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
object, each time calling
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
with that model's logit column as a numeric `distance` Every matching
option is passed explicitly, so a future change of a MatchIt default
cannot alter the results. All models are matched, whether or not they
passed the log-loss gate; the gate is applied later by the wrapper
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md).

## Usage

``` r
match_each_model(
  model_set,
  data,
  replace = TRUE,
  ratio = 1,
  m_order = "largest",
  caliper = NULL,
  oracle_logit = NULL,
  bootstrap_samples = 0,
  cores = 1,
  seed = 1,
  common_random_order = FALSE,
  clip_bootstrap = FALSE
)

# S3 method for class 'match_set'
print(x, ...)
```

## Arguments

- model_set:

  A `propensity_model_set` from
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md).

- data:

  The data frame the model set was fitted on (same rows, same order).
  The matches are stored with the row names of `data` as ids; the data
  itself is not stored.

- replace:

  `TRUE` (matching with replacement; the default of this function) or
  `FALSE` (without replacement; the default of
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  with a formula). Both are supported; see
  [`vignette("definitions", package = "propsens")`](https://mostafa-abbas.github.io/propsens/articles/definitions.md)
  for what each targets.

- ratio:

  Number of untreated partners per treated patient. Only `1` is
  supported in this version.

- m_order:

  Order in which treated patients are matched; one of MatchIt's
  `"largest"`, `"smallest"`, `"data"`, `"random"`, `"closest"`,
  `"farthest"`.

- caliper:

  `NULL` (default, no caliper) or a positive number: the caliper in
  standard deviations of the logit.

- oracle_logit:

  Internal, for simulation studies only: a numeric vector of length n
  with the true logit, for example the `true_logit` column of
  [`simulate_example_data()`](https://mostafa-abbas.github.io/propsens/reference/simulate_example_data.md).
  When given, one extra entry `oracle` is matched with the same options.
  Not needed in an analysis.

- bootstrap_samples:

  Number of bootstrap samples `B` for the sampling-instability block
  (default 0, none).

- cores:

  Number of cores for the matching tasks (default 1).

- seed:

  Integer seed from which every other seed is derived.

- common_random_order:

  Internal, for simulation studies only, and only for
  `m_order = "random"`: `FALSE` (default) gives each model its own
  random order; `TRUE` gives every model (and the oracle) the same
  random order, so that the models differ only in their distances, not
  in the order in which treated patients are matched.

- clip_bootstrap:

  `FALSE` (default) or `TRUE`: clip the bootstrap refits' scores to
  `[clip, 1 - clip]` with the model set's `clip`, as the cross-fitted
  scores are. See the section "Bootstrap block".

- x:

  A `match_set` object.

- ...:

  Not used.

## Value

An object of class `"match_set"`: a list with

- `matches`: one data frame per model (`treated_id`, `control_id`,
  `weight`, `pair_id`), named by model id;

- `weights`: n x M matrix of MatchIt weights, rows named by the row ids;

- `treated_kept`: named vector, share of treated patients that were
  matched;

- `calls`: named character vector with the deparsed `matchit()` call of
  each model; the call names the data and the model set as the caller
  wrote them and re-evaluates to the same matches;

- `oracle`: `NULL`, or a list with `matches`, `weights`, `treated_kept`
  and `call`;

- `bootstrap`: `NULL`, or a list with `rank`, `pairs`, `present` (n x B
  matrices) and `skipped` (data frame of skipped samples with the
  reason) and `design_rule_samples` (the samples whose refit used the
  design-matrix rule for a level absent from the sample) and
  `prediction_fallback_samples` (the samples whose refit used it because
  prediction failed);

- `models`: the model table of the model set (`model_id`, `algorithm`,
  `seed_index`, `passed_logloss_gate`), used by the diagnostics to apply
  the quality gate;

- `model_ids`, `treatment`, `row_ids`, `settings` and `call`.

## Matching options

The call is
`MatchIt::matchit(formula, data, method = "nearest", estimand = "ATT", distance = <logit column>, replace = replace, ratio = 1, m.order = m_order, caliper = caliper, std.caliper = TRUE, normalize = FALSE)`.
The order `"largest"` is MatchIt's default for the ATT when a numeric
distance is supplied: the treated patient with the largest score is
matched first. A caliper is given in standard deviations of the logit
(`std.caliper = TRUE`); it drops treated patients, so the share of
treated kept is reported.

## Weights

`normalize = FALSE` is passed so that MatchIt's weights are the raw
matching weights: every matched treated patient has weight 1, every
matched untreated patient has weight equal to the number of pairs it
belongs to, and unmatched patients have weight 0. The sum of the weights
of the untreated patients is therefore the number of matched treated
patients, also with replacement. (MatchIt's default `normalize = TRUE`
rescales those weights to a mean of 1, which changes no weighted mean or
balance statistic but would make that sum equal the number of distinct
untreated partners.) Each row of a match table carries
`weight = 1 / ratio`, which is 1 in this version.

## Bootstrap block

With `bootstrap_samples = B > 0` the function draws `B` nonparametric
row bootstrap samples of the data (one derived seed per sample). A
sample with fewer than 10 treated or 10 untreated rows is skipped,
listed in `bootstrap$skipped` and counted in a message. For each sample,
logistic regression with the model set's formula (as written, with any
interactions or transformations) is fitted in-sample, without
cross-fitting: this block measures the sampling variability of the
fitted model, not its out-of-sample accuracy (see
[`bootstrap_rank_instability()`](https://mostafa-abbas.github.io/propsens/reference/bootstrap_rank_instability.md)).
The refitted model predicts logits for every ORIGINAL patient; their
percentile ranks form `bootstrap$rank` (n x B). These logits are not
clipped by default, while the cross-fitted scores are clipped to
`[clip, 1 - clip]`; under near-separation many cross-fitted scores can
then share the clipped value (tied ranks) while their bootstrap ranks
stay distinct. `clip_bootstrap = TRUE` clips the bootstrap scores in the
same way. The default (`FALSE`) keeps the method as it was studied. When
a factor level of the original data is absent from the sample, the refit
uses the full-data design matrix and scores those patients as the
reference level (the rule of the cross-fitted logistic model); the
samples concerned are listed in `bootstrap$design_rule_samples`. When
every level is present but prediction from the refitted model fails (a
rank-deficient fit that did not converge), the refit uses the same
design-matrix rule; those samples are listed in
`bootstrap$prediction_fallback_samples`. No sample is skipped for this
reason. The sample is then matched with the same options, rows are
mapped back to original ids, and for every original treated patient
present in the sample the untreated partner matched to its first
occurrence is stored in `bootstrap$pairs` (n x B, `NA` when absent or
unmatched). `bootstrap$present` (n x B, logical) says whether each
original patient is in each sample. A skipped sample keeps its column,
filled with `NA`, so that every matrix is n x B.

## Seeds and parallel runs

Per-task seeds are derived in the main process from `seed`, passed to
each worker and applied with
[`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html);
`m_order = "random"` and the bootstrap draws are therefore reproducible
and identical for `cores = 1` and `cores > 1`. The matching tasks, the
bootstrap draws and the bootstrap matching each have their own family of
derived seeds, separate from the model fits; with the default order
`"largest"` matching uses no random numbers. Parallel runs use
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
Unix; on Windows the run is serial.

## Examples

``` r
data("lalonde", package = "MatchIt")
model_set <- fit_propensity_models(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, seeds = 1, n_trees = 100,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
)
with_replacement <- match_each_model(model_set, lalonde, replace = TRUE)
with_replacement
#> Match set: 3 model(s), 614 patients, 185 treated
#> Matching: 1:1 nearest neighbour on the logit, ATT, with replacement, order "largest", caliper none
#>               model_id treated_kept distinct_controls
#>  logistic_regression_1            1                85
#>        random_forest_1            1                77
#>    gradient_boosting_1            1                69
head(with_replacement$matches$logistic_regression_1)
#>   treated_id control_id weight pair_id
#> 1       NSW1    PSID157      1       1
#> 2       NSW2    PSID224      1       2
#> 3       NSW3    PSID335      1       3
#> 4       NSW4     PSID15      1       4
#> 5       NSW5    PSID401      1       5
#> 6       NSW6    PSID134      1       6
without_replacement <- match_each_model(model_set, lalonde, replace = FALSE,
                                        bootstrap_samples = 5)
without_replacement
#> Match set: 3 model(s), 614 patients, 185 treated
#> Matching: 1:1 nearest neighbour on the logit, ATT, without replacement, order "largest", caliper none
#>               model_id treated_kept distinct_controls
#>  logistic_regression_1            1               185
#>        random_forest_1            1               185
#>    gradient_boosting_1            1               185
#> Bootstrap: 5 logistic refits (0 sample(s) skipped)
dim(without_replacement$bootstrap$rank)
#> [1] 614   5
```
