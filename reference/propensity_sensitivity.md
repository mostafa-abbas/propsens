# Propensity model sensitivity of a matched design

The one-call diagnostic. It refits the propensity model with every
algorithm of the model set, re-matches with the analyst's matching
options, applies the quality gate, and reports rank disagreement (the
primary summary), pair instability (the second summary of the same gated
set: 1 minus the pair recurrence across models; exploratory), the
bootstrap rank instability of the logistic model, control spread and
pair recurrence, and study-level balance spread. No outcome is used.
Rank disagreement is a study-level diagnostic; read it for the study as
a whole, not per patient. The diagnostics are evaluated in a manuscript
in preparation. When fewer than two models pass the quality gate,
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) say so first: the
diagnostic is then outside the setting where it was evaluated.

## Usage

``` r
propensity_sensitivity(
  x,
  data = NULL,
  algorithms = c("logistic_regression", "bart", "grf_probability"),
  seeds = 1:3,
  folds = 5,
  cross_fit = TRUE,
  n_trees = NULL,
  quality_tolerance = 0.05,
  balance_rule = c("relative", "absolute", "relative_only"),
  balance_margin = 0.1,
  balance_cap = 0.25,
  balance_threshold = 0.1,
  bootstrap_samples = 20,
  covariates = NULL,
  replace = FALSE,
  m_order = "largest",
  caliper = NULL,
  cores = 1,
  seed = 1
)

# S3 method for class 'propensity_sensitivity'
print(x, digits = 3, ...)

# S3 method for class 'propensity_sensitivity'
summary(object, ...)

# S3 method for class 'summary.propensity_sensitivity'
print(x, digits = 3, ...)

# S3 method for class 'propensity_sensitivity'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `matchit` object (from
  [`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html))
  or a two-sided formula such as `treat ~ x1 + x2`. Only baseline
  covariates belong in the formula; exclude outcomes, identifiers and
  post-treatment variables (for example, `treat ~ . - death - id`).

- data:

  The data frame. Required with a formula; with a `matchit` object it is
  looked up from the object's call when `NULL`.

- algorithms:

  The propensity model set. The default is the set evaluated in a
  manuscript in preparation: logistic regression, BART and the GRF
  probability forest (with the default `seeds`, three seeds for each
  tree algorithm: seven models; packages 'dbarts' and 'grf'). A faster
  set is available with
  `algorithms = c("logistic_regression", "random_forest", "gradient_boosting")`;
  it is not part of that evaluation. See
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
  for every option and the fixed settings.

- seeds, folds, cross_fit, n_trees, quality_tolerance:

  Passed to
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md).
  `n_trees` applies to gradient boosting only.

- balance_rule:

  The balance part of the quality gate: `"relative"` (default),
  `"relative_only"` (the relative rule without the cap) or `"absolute"`.
  See Steps.

- balance_margin, balance_cap:

  Relative rule: a model fails if its post-match maximum absolute
  standardized mean difference on main terms is more than
  `balance_margin` above the best model's (default 0.10) or above
  `balance_cap` (default 0.25; not used by `"relative_only"`).

- balance_threshold:

  Absolute rule: largest allowed post-match maximum absolute
  standardized mean difference on main terms (default 0.10). Used only
  with `balance_rule = "absolute"`.

- bootstrap_samples:

  Number of logistic bootstrap refits, for the bootstrap rank
  instability and the bootstrap pair recurrence (default 20; 0 for
  none).

- covariates:

  Data frame of covariates for the control spread and the balance
  statistics. Its rows must be the same patients in the same order as
  `data`; when both have row names other than the default and they
  differ, the function stops. `NULL` (default) uses the source columns
  of the formula's right-hand side, from `data`, resolved with R's
  formula rules: `.` means every other column and exclusions such as
  `- death` are respected; `I(age^2)` gives `age`, `age:diabetes` gives
  `age` and `diabetes`.

- replace, m_order, caliper:

  Matching options for the formula interface (without replacement, order
  `"largest"`, no caliper by default; the caliper is in standard
  deviations of each model's logit). With a `matchit` object they are
  read from the object and must not be given.

- cores:

  Number of cores (default 1).

- seed:

  Integer seed from which every other seed is derived.

- digits:

  Significant digits for printed numbers.

- ...:

  Not used.

- object:

  A `propensity_sensitivity` object.

- row.names, optional:

  Ignored; present for the generic.

## Value

An object of class `"propensity_sensitivity"`: a list with `patients`
(one row per patient: `id`, `treated`, `rank_disagreement`,
`rank_disagreement_between`, `rank_disagreement_within`,
`logit_disagreement`, `control_spread`, `control_spread_scaled`,
`pair_recurrence_models`, `pair_recurrence_bootstrap`,
`outside_support`, `logistic_score`; the treated-only measures are `NA`
for untreated patients), `models` (one row per model plus the reference
row), `summary` (a table of the study-level numbers with and without the
gate, plus flags), `quality_gate_left_fewer_than_two_models`,
`balance_spread` (spread tables with and without the gate),
`bootstrap_rank_instability` (see
[`bootstrap_rank_instability()`](https://mostafa-abbas.github.io/propsens/reference/bootstrap_rank_instability.md);
`NULL` without bootstrap refits; its study-level value is also the row
`bootstrap_rank_instability_treated` of the summary table; the summary
table also has `pair_instability_models_mean_treated`, 1 minus
`pair_recurrence_models_mean_treated`), `settings` and `components` (the
model set and the match set). `components$model_set` is the fitted
[`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
object, with every model's cross-fitted `probability`, `logit` and
`rank`, so that its scores can be reused (for example in another
estimator) without refitting; `summary$models_used` names the models in
use.

## Steps

1.  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
    with the formula (cross-fitted scores by default).

2.  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md)
    with the matching options and the logistic bootstrap block.

3.  [`balance_spread()`](https://mostafa-abbas.github.io/propsens/reference/balance_spread.md)
    for every model.

4.  Final quality gate:
    `passed_quality_gate = passed_logloss_gate & smd_max_main <= balance_limit`.
    With the default `balance_rule = "relative"` the limit is the
    smallest max \|SMD\| on main terms among the refitted models plus
    `balance_margin` (0.10), and never above `balance_cap` (0.25): a
    model is excluded only when its balance is clearly worse than the
    best model's. With `balance_rule = "relative_only"` the limit is the
    best model's value plus `balance_margin`, without the cap; use it as
    a sensitivity analysis when the cap excludes every model (many
    covariates). With `balance_rule = "absolute"` the limit is
    `balance_threshold` (0.10). If fewer than two models pass, every
    model is used and `quality_gate_left_fewer_than_two_models` is
    `TRUE` (applied to the final set).

5.  [`rank_disagreement()`](https://mostafa-abbas.github.io/propsens/reference/rank_disagreement.md),
    [`control_spread()`](https://mostafa-abbas.github.io/propsens/reference/control_spread.md),
    [`pair_recurrence()`](https://mostafa-abbas.github.io/propsens/reference/pair_recurrence.md)
    and the balance spread on the final set (`summary` column
    `with_gate`) and on every model (`without_gate`).

## Reading a matchit object

From a `matchit` object the function reads the formula, the method,
`ratio`, `replace`, `m.order` (from the call; MatchIt's default for the
ATT with a propensity score is `"largest"`) and the caliper. Only 1:1
nearest-neighbour matching on a propensity score for the ATT is
supported, with or without replacement and with or without a caliper on
the score; anything else stops with a clear error. MatchIt applies a
caliper in standard deviations of the analyst's distance (over all
patients); the same number of standard deviations is applied to the
logit of every model here. When `data` is `NULL` the data is found by
evaluating the object's `data` argument in the calling frame.

## Reference model

The analyst's own score (`x$distance`) is kept as the reference model,
an extra row of `models` with `model_id = "reference"`. Its ranks come
from `x$distance` (ranks do not depend on the scale) and its balance
statistics from `x$weights`. Its log-loss, Brier score and c-statistic
are `NA`: it is an in-sample score from the user's object and is not
comparable with the cross-fitted scores. Its quality gate is judged on
balance only; when it fails, the row is flagged "reference model outside
the equally-reasonable set". It is never one of the models in use.
`summary` reports its mean absolute rank difference from the models in
use, over treated patients: the mean over treated `i` and models `m` of
`|r_ref(i) - r_m(i)|`.

## Examples

``` r
data("lalonde", package = "MatchIt")
matched <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                            data = lalonde, distance = "glm", link = "linear.logit",
                            replace = TRUE)
# A fast model set for this example; the default set is slower (see below).
sensitivity <- propensity_sensitivity(
  matched, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
  seeds = 1, n_trees = 100, bootstrap_samples = 5
)
sensitivity
#> Propensity model sensitivity
#> WARNING: fewer than two models passed the quality gate; the diagnostic is outside the setting where it was evaluated.
#> Every model is used below, so the values do not describe a set of equally good models. Do not read them as a validated result.
#> 1. Quality gate: 3 of 3 models in use (3 algorithms: logistic regression, random forest, gradient boosting); log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
#>    Fewer than two models passed, so every model is used. Failed:
#>      logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      gradient_boosting_1: max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#> 2. Rank disagreement (study mean over treated): 0.0687; median between/within ratio NA (exploratory; near 1 = seed noise)
#>    Pair instability across the models in use (study mean over treated; exploratory): 0.6
#>    Bootstrap rank instability of the logistic model (study mean over treated): 0.0361
#> 3. Balance spread (range of max |SMD| main terms across the models in use): 0.125 (0.16 to 0.285)
#> Mean scaled control spread over treated (descriptive): 1.19 (1 = controls differ as much as a typical match differs from its treated patient)
#> Treated kept: 100 % (lowest model in use)
#> Reference model (glm, linear.logit): passes the balance gate; mean absolute rank difference from the models in use: 0.0766
summary(sensitivity)
#> WARNING: fewer than two models passed the quality gate; the diagnostic is outside the setting where it was evaluated.
#> Every model is used below, so the values do not describe a set of equally good models. Do not read them as a validated result.
#> 1. Quality gate: 3 of 3 models in use (3 algorithms: logistic regression, random forest, gradient boosting); log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
#>    Fewer than two models passed, so every model is used. Failed:
#>      logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      gradient_boosting_1: max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#> 2. Rank disagreement (study mean over treated): 0.0687; median between/within ratio NA (exploratory; near 1 = seed noise)
#>    Pair instability across the models in use (study mean over treated; exploratory): 0.6
#>    Bootstrap rank instability of the logistic model (study mean over treated): 0.0361
#> 3. Balance spread (range of max |SMD| main terms across the models in use): 0.125 (0.16 to 0.285)
#> 
#> Study-level numbers (with the quality gate / without it):
#>                                         quantity with_gate without_gate
#>                            number_of_models_used    3.0000       3.0000
#>                   mean_rank_disagreement_treated    0.0687       0.0687
#>                       mean_rank_disagreement_all    0.0768       0.0768
#>              median_between_within_ratio_treated        NA           NA
#>                  mean_logit_disagreement_treated    0.7330       0.7330
#>               mean_control_spread_scaled_treated    1.1900       1.1900
#>                p90_control_spread_scaled_treated    1.8400       1.8400
#>                               smd_max_main_range    0.1250       0.1250
#>                                  smd_max_main_sd    0.0686       0.0686
#>                           smd_max_extended_range    0.2150       0.2150
#>                              smd_max_extended_sd    0.1080       0.1080
#>                            energy_distance_range    0.0398       0.0398
#>                               energy_distance_sd    0.0213       0.0213
#>                                     ks_max_range    0.2220       0.2220
#>                                        ks_max_sd    0.1260       0.1260
#>              pair_recurrence_models_mean_treated    0.4000       0.4000
#>             pair_instability_models_mean_treated    0.6000       0.6000
#>           pair_recurrence_bootstrap_mean_treated    0.6200       0.6200
#>                               treated_kept_share    1.0000       1.0000
#>                           treated_kept_share_min    1.0000       1.0000
#>               bootstrap_rank_instability_treated    0.0361       0.0361
#>  reference_mean_absolute_rank_difference_treated    0.0766       0.0766
#> 
#> Models:
#>               model_id                     algorithm logloss smd_max_main
#>  logistic_regression_1           logistic_regression   0.412        0.271
#>        random_forest_1                 random_forest   0.363        0.160
#>    gradient_boosting_1             gradient_boosting   0.351        0.285
#>              reference reference (glm, linear.logit)      NA        0.224
#>  smd_max_extended energy_distance ks_max treated_kept passed_quality_gate
#>             0.459          0.0953  0.357            1               FALSE
#>             0.244          0.0555  0.141            1                TRUE
#>             0.377          0.0885  0.135            1               FALSE
#>             0.233          0.0806  0.335            1                TRUE
#>  in_use
#>    TRUE
#>    TRUE
#>    TRUE
#>   FALSE
#>                                                                                                                                                                            gate_note
#>  log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap; used anyway: fewer than two models passed
#>                                                                                                                                                                               passed
#>                                                 max |SMD| main terms more than balance_margin above the best model's or above balance_cap; used anyway: fewer than two models passed
#>                                                                  in-sample score from the user's object, not comparable (log-loss not computed); gate judged on balance only: passed
#> 
#> Quality gate: log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
#> Matching: 1:1 nearest neighbour on the logit, ATT, with replacement, order "largest", caliper none
head(as.data.frame(sensitivity))
#>     id treated rank_disagreement rank_disagreement_between
#> 1 NSW1       1        0.04782657                0.04782657
#> 2 NSW2       1        0.12597696                0.12597696
#> 3 NSW3       1        0.11200054                0.11200054
#> 4 NSW4       1        0.01447589                0.01447589
#> 5 NSW5       1        0.05719102                0.05719102
#> 6 NSW6       1        0.06427725                0.06427725
#>   rank_disagreement_within logit_disagreement control_spread
#> 1                       NA         0.06197223      3.5768753
#> 2                       NA         1.07546033      4.8896518
#> 3                       NA         1.25276825      1.7243229
#> 4                       NA         1.94151550      3.8813720
#> 5                       NA         0.84889473      0.9260135
#> 6                       NA         1.26714331      2.1697726
#>   control_spread_scaled pair_recurrence_models pair_recurrence_bootstrap
#> 1             1.3842876              0.3333333                 0.3333333
#> 2             1.8923457              0.3333333                 1.0000000
#> 3             0.6673307              0.6666667                 0.2500000
#> 4             1.5021310              0.3333333                 1.0000000
#> 5             0.3583768              0.3333333                 0.3333333
#> 6             0.8397244              0.3333333                 0.2500000
#>   outside_support logistic_score
#> 1           FALSE      0.3950590
#> 2           FALSE     -1.3019135
#> 3           FALSE      0.5733472
#> 4           FALSE      1.2414886
#> 5           FALSE      0.8807265
#> 6           FALSE      0.9013876
# \donttest{
# The defaults: logistic regression, BART x 3, GRF probability forest x 3 (about a minute).
propensity_sensitivity(matched)
#> Propensity model sensitivity
#> WARNING: fewer than two models passed the quality gate; the diagnostic is outside the setting where it was evaluated.
#> Every model is used below, so the values do not describe a set of equally good models. Do not read them as a validated result.
#> 1. Quality gate: 7 of 7 models in use (3 algorithms: logistic regression, BART, GRF); log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
#>    Fewer than two models passed, so every model is used. Failed:
#>      logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      bart_1: max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      bart_2: max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      bart_3: max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      grf_probability_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      grf_probability_2: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#>      grf_probability_3: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
#> 2. Rank disagreement (study mean over treated): 0.0493; median between/within ratio 10.1 (exploratory; near 1 = seed noise)
#>    Pair instability across the models in use (study mean over treated; exploratory): 0.659
#>    Bootstrap rank instability of the logistic model (study mean over treated): 0.0462
#> 3. Balance spread (range of max |SMD| main terms across the models in use): 0.297 (0.262 to 0.559)
#> Mean scaled control spread over treated (descriptive): 1.24 (1 = controls differ as much as a typical match differs from its treated patient)
#> Treated kept: 100 % (lowest model in use)
#> Reference model (glm, linear.logit): passes the balance gate; mean absolute rank difference from the models in use: 0.0903
# }
```
