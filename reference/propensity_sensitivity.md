# Propensity model sensitivity of a matched design

How much does a matched design depend on which of several similarly good
propensity models is used? The function refits the propensity model with
every algorithm of the model set, repeats the analyst's matching with
each model, keeps the models that pass a quality gate (the kept models),
and describes how much the kept models disagree. No outcome is used.

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
  seed = 1,
  min_treated_kept = NULL,
  clip_bootstrap = FALSE
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
  post-treatment variables (for example, `treat ~ . - death - id`). The
  treatment may be numeric 0/1, logical or a two-level factor (see
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)).

- data:

  The data frame. Required with a formula; with a `matchit` object it is
  looked up from the object's call when `NULL`.

- algorithms:

  The propensity model set. The default is logistic regression, BART and
  the GRF probability forest (with the default `seeds`, three seeds for
  each tree algorithm: seven models; packages 'dbarts' and 'grf'). A
  faster set is available with
  `algorithms = c("logistic_regression", "random_forest", "gradient_boosting")`
  (packages 'randomForest' and 'gbm'); it has been studied less. See
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
  for every option and the fixed settings.

- seeds:

  `1:k`, the number of seeds per tree algorithm (default `1:3`). Only
  the number is used; set the randomness with `seed`.

- folds, cross_fit, n_trees, quality_tolerance:

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
  instability and the bootstrap pair recurrence: 0 for none, otherwise
  at least 2 (default 20).

- covariates:

  Data frame of covariates for the control spread and the balance
  statistics. Its rows must be the same patients in the same order as
  `data`. The row names are checked when `covariates` has row names
  other than the default: if `data` has them too, they must be
  identical; if `data` has the default row names (1 to n), those of
  `covariates` must be "1" to "n" in order (a reordered copy of `data`
  is an error), and other row names are an error. When `covariates` has
  the default row names (for example a data frame built column by
  column) the order cannot be checked and is your responsibility. `NULL`
  (default) uses the source columns of the formula's right-hand side,
  from `data`, resolved with R's formula rules: `.` means every other
  column and exclusions such as `- death` are respected; `I(age^2)`
  gives `age`, `age:diabetes` gives `age` and `diabetes`. Character
  columns are an error (convert them to factors).

- replace, m_order, caliper:

  Matching options for the formula interface (without replacement, order
  `"largest"`, no caliper by default; the caliper is in standard
  deviations of each model's logit). With a `matchit` object they are
  read from the object and must not be given.

- cores:

  Number of cores (default 1). Run time grows with the number of
  patients and covariates (about a minute for 614 patients and 7
  covariates; about 16 minutes for 5,735 patients and 50 covariates on
  one core); on Unix, `cores` runs the models in parallel.

- seed:

  Integer seed from which every other seed is derived (fold split,
  models, bootstrap). See the section "Seeds".

- min_treated_kept:

  `NULL` (default, off) or a share between 0 and 1: with a caliper, a
  model that keeps a smaller share of the treated patients fails the
  quality gate and is left out before the relative balance limit is
  computed.

- clip_bootstrap:

  `FALSE` (default) or `TRUE`: clip the bootstrap refits' scores to
  `[0.001, 0.999]` as the cross-fitted scores are. The default keeps the
  method as it was studied; see
  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md).

- digits:

  Significant digits for printed numbers.

- ...:

  Not used.

- object:

  A `propensity_sensitivity` object.

- row.names, optional:

  Ignored; present for the generic.

## Value

An object of class `"propensity_sensitivity"`: a list with

- `patients`: one row per patient (`id`, `treated`, `rank_disagreement`,
  `rank_disagreement_between`, `rank_disagreement_within`,
  `logit_disagreement`, `control_spread`, `control_spread_scaled`,
  `pair_recurrence_models`, `pair_recurrence_bootstrap`,
  `outside_support`, `logistic_score`; the treated-only measures are
  `NA` for untreated patients). These columns describe the spread over
  patients; they are not per-patient flags.

- `models`: one row per model (fit, balance, `treated_kept`, the gate
  columns, `in_use` = kept, `gate_note`), plus the row of your own model
  for a `matchit` object.

- `summary`: a list with `table` (the study-level numbers with the gate,
  `with_gate`, and without it, `without_gate`),
  `quality_gate_left_fewer_than_two_models`, `reference_model_not_kept`,
  `models_used` (the kept models), `kept_algorithms`,
  `kept_algorithms_count`, `kept_models_one_algorithm`,
  `pair_recurrence_floor_models`, `pair_instability_maximum` and
  `pair_recurrence_expected_ceiling_bootstrap_uncorrected`.

- `quality_gate_left_fewer_than_two_models` and `kept_algorithms_count`
  (also at the top level, for scripts).

- `balance_spread`: spread tables with and without the gate.

- `bootstrap_rank_instability`: see
  [`bootstrap_rank_instability()`](https://mostafa-abbas.github.io/propsens/reference/bootstrap_rank_instability.md);
  `NULL` without bootstrap refits.

- `settings` and `components` (`model_set`, the fitted
  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
  object with every model's cross-fitted `probability`, `logit` and
  `rank`, so that the scores can be reused without refitting; and
  `match_set`).

## Details

Every summary is study-level: read it for the study as a whole. Rank
disagreement is the primary summary. Pair instability, rank disagreement
with equal weight per algorithm and the between/within-algorithm ratio
are exploratory. Bootstrap rank instability, balance spread and control
spread are descriptive. The per-patient columns of `patients` (and of
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)) describe
how the disagreement is spread over the patients; they are not flags for
individual patients. No thresholds for a "large" value exist yet.

These diagnostics have been studied in the author's simulations and
benchmark analyses; a methods paper is in preparation. No peer-reviewed
evaluation exists yet.

Two situations raise an R warning and a block in
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html): fewer than two
models pass the quality gate (every model is then used, so there is no
set of similarly good models to compare), and the kept models are seeds
of a single algorithm (rank disagreement and pair instability then
measure seed noise within that algorithm).

## Steps

1.  All arguments, the treatment and the covariates are checked before
    any model is fitted.

2.  [`fit_propensity_models()`](https://mostafa-abbas.github.io/propsens/reference/fit_propensity_models.md)
    with the formula (cross-fitted scores by default).

3.  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md)
    with the matching options and the logistic bootstrap block.

4.  [`balance_spread()`](https://mostafa-abbas.github.io/propsens/reference/balance_spread.md)
    for every model.

5.  Quality gate:
    `passed_quality_gate = passed_logloss_gate & smd_max_main <= balance_limit`.
    With the default `balance_rule = "relative"` the limit is the
    smallest max \|SMD\| on main terms among the refitted models plus
    `balance_margin` (0.10), and never above `balance_cap` (0.25): a
    model is excluded only when its balance is clearly worse than the
    best model's. With `balance_rule = "relative_only"` the limit is the
    best model's value plus `balance_margin`, without the cap; use it as
    a sensitivity analysis when the cap excludes every model (many
    covariates). With `balance_rule = "absolute"` the limit is
    `balance_threshold` (0.10). With `min_treated_kept` set (off by
    default), a model that keeps a smaller share of the treated patients
    fails the gate and does not enter the best-model value of the
    relative rule. If fewer than two models pass, every model is used
    and `quality_gate_left_fewer_than_two_models` is `TRUE`.

6.  [`rank_disagreement()`](https://mostafa-abbas.github.io/propsens/reference/rank_disagreement.md),
    [`control_spread()`](https://mostafa-abbas.github.io/propsens/reference/control_spread.md),
    [`pair_recurrence()`](https://mostafa-abbas.github.io/propsens/reference/pair_recurrence.md)
    and the balance spread on the kept models (`summary` column
    `with_gate`) and on every model (`without_gate`).

## Seeds

The fold split and every tree model depend on `seed`. On small data the
kept models and every summary can change noticeably with `seed` (on
lalonde, rank disagreement varied several-fold between seeds). Rerun
with two or three values of `seed` and report all of them;
[`seed_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/seed_sensitivity.md)
does this. Rank disagreement is a standard deviation over the kept
models, so it also depends on how many seeds each kept algorithm
contributes; "rank disagreement, equal weight per algorithm" counts each
algorithm once.

## What the summaries compare under a caliper

With a caliper each model keeps a different subset of the treated
patients, shown by [`print()`](https://rdrr.io/r/base/print.html) for
each kept model. The standardized mean differences compare the matched
treated patients with their matched untreated partners; the energy
distance (`cobalt::bal.init(estimand = "ATT")`) also measures the
distance to the whole treated group, so it also grows when treated
patients are dropped. Rank disagreement averages over all treated
patients, matched or not; pair instability and control spread use the
treated patients matched by at least one kept model. The gate does not
look at the share of treated kept unless `min_treated_kept` is set: a
model that keeps few treated patients can reach an easy balance and then
set the relative limit for the others.

## Estimand

Matching is 1:1 nearest-neighbour on the logit of each model's score,
for the average treatment effect in the treated (ATT). Without
replacement (the default of the formula interface and the more common
applied choice) the result depends on the matching order, and late
matches can be poor when untreated patients are not much more numerous
than treated ones, so the matched estimate can keep some bias for the
ATT. With replacement every treated patient gets its nearest untreated
patient; under overlap and a correctly specified propensity model this
targets the ATT, at the price of reused untreated patients and a
variance that must account for it. With a caliper, treated patients
without an eligible partner are dropped and the estimand is the effect
in the retained treated patients, not the ATT.

## Reading a matchit object

From a `matchit` object the function reads the formula, the method,
`ratio`, `replace`, `m.order` (from the call; MatchIt's default for the
ATT with a propensity score is `"largest"`) and the caliper. Only 1:1
nearest-neighbour matching on a propensity score for the ATT is
supported, with or without replacement and with or without a caliper on
the score; anything else (including `exact`, `mahvars`, `antiexact`,
`s.weights`, `reuse.max`, `unit.id` and calipers on covariates) stops
with a clear error. MatchIt applies a caliper in standard deviations of
the analyst's distance (over all patients); the same number of standard
deviations is applied to the logit of every model here. When `data` is
`NULL` the data is found by evaluating the object's `data` argument in
the calling frame.

## Your own model

The analyst's own score (`x$distance`) is kept as an extra row of
`models` with `model_id = "reference"`. Its ranks come from `x$distance`
(ranks do not depend on the scale) and its balance statistics from
`x$weights`. Its log-loss, Brier score and c-statistic are `NA`: it is
an in-sample score from the user's object and is not comparable with the
cross-fitted scores. It is judged on balance only, as "within the
balance limit" or "above the balance limit"; it is never one of the
models compared, so it is never "kept". `summary` reports its mean
absolute rank difference from the kept models, over treated patients:
the mean over treated `i` and models `m` of `|r_ref(i) - r_m(i)|`.

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
#> Warning: Fewer than two models passed the quality gate, so there is no set of similarly good models to compare. Every model is used below; read the values as a description of all models, not of a kept set.
sensitivity
#> Propensity model sensitivity
#> WARNING: fewer than two models passed the quality gate, so there is no set of
#> similarly good models to compare. Every model is used below; read the values
#> as a description of all models, not of a kept set.
#> 1. Quality gate: 3 of 3 models in use (3 algorithms: logistic regression,
#>    random forest, gradient boosting); log-loss within 5 % of the best;
#>    relative balance rule: max |SMD| main terms at most 0.10 above the best
#>    model's and at most 0.25 (limit here 0.25).
#>    Fewer than two models passed, so every model is used. Failed:
#>      logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
#>        |SMD| main terms more than balance_margin above the best model's or
#>        above balance_cap
#>      gradient_boosting_1: max |SMD| main terms more than balance_margin above
#>        the best model's or above balance_cap
#>    Treated patients kept by each model: logistic_regression_1 100 %,
#>      random_forest_1 100 %, gradient_boosting_1 100 %
#> 2. Rank disagreement (primary; study mean over treated): 0.0687
#>    Equal weight per algorithm (exploratory): 0.0687
#>    Between/within-algorithm ratio (exploratory; median over treated): NA
#>    Pair instability (exploratory; study mean over treated): 0.6 of a possible
#>      0.667 with 3 models. It is high whenever many untreated patients have
#>      similar scores, even when the models agree.
#>    Bootstrap rank instability of the logistic model (descriptive): 0.0361
#> 3. Balance spread (descriptive; range of max |SMD| main terms over the models
#>    in use): 0.125 (0.16 to 0.285)
#> Control spread (descriptive; mean scaled over treated): 1.19; 1 = the partners
#>    chosen by different models differ as much as a typical matched pair (scale:
#>    the matches of logistic_regression_1)
#> Your own model (glm, linear.logit): within the balance limit (judged on
#>    balance only; it is not one of the compared models); mean absolute rank
#>    difference from the models in use: 0.0766
#> Matching: 1:1 nearest neighbour on the logit, with replacement, order
#>    "largest", caliper none; target: ATT.
summary(sensitivity)
#> Propensity model sensitivity: summary
#> WARNING: fewer than two models passed the quality gate, so there is no set of
#> similarly good models to compare. Every model is used below; read the values
#> as a description of all models, not of a kept set.
#> 1. Quality gate: 3 of 3 models in use (3 algorithms: logistic regression,
#>    random forest, gradient boosting); log-loss within 5 % of the best;
#>    relative balance rule: max |SMD| main terms at most 0.10 above the best
#>    model's and at most 0.25 (limit here 0.25).
#>    Fewer than two models passed, so every model is used. Failed:
#>      logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
#>        |SMD| main terms more than balance_margin above the best model's or
#>        above balance_cap
#>      gradient_boosting_1: max |SMD| main terms more than balance_margin above
#>        the best model's or above balance_cap
#>    Treated patients kept by each model: logistic_regression_1 100 %,
#>      random_forest_1 100 %, gradient_boosting_1 100 %
#> 2. Rank disagreement (primary; study mean over treated): 0.0687
#>    Equal weight per algorithm (exploratory): 0.0687
#>    Between/within-algorithm ratio (exploratory; median over treated): NA
#>    Pair instability (exploratory; study mean over treated): 0.6 of a possible
#>      0.667 with 3 models. It is high whenever many untreated patients have
#>      similar scores, even when the models agree.
#>    Bootstrap rank instability of the logistic model (descriptive): 0.0361
#> 3. Balance spread (descriptive; range of max |SMD| main terms over the models
#>    in use): 0.125 (0.16 to 0.285)
#> Control spread (descriptive; mean scaled over treated): 1.19; 1 = the partners
#>    chosen by different models differ as much as a typical matched pair (scale:
#>    the matches of logistic_regression_1)
#> Your own model (glm, linear.logit): within the balance limit (judged on
#>    balance only; it is not one of the compared models); mean absolute rank
#>    difference from the models in use: 0.0766
#> Matching: 1:1 nearest neighbour on the logit, with replacement, order
#>    "largest", caliper none; target: ATT.
#> 
#> Study-level numbers (every model is used, so one column):
#>  quantity                                        label       value 
#>  number_of_models_used                           count       3     
#>  mean_rank_disagreement_treated                  primary     0.0687
#>  mean_rank_disagreement_equal_weight_treated     exploratory 0.0687
#>  mean_rank_disagreement_all                      descriptive 0.0768
#>  median_between_within_ratio_treated             exploratory NA    
#>  mean_logit_disagreement_treated                 descriptive 0.733 
#>  mean_control_spread_scaled_treated              descriptive 1.19  
#>  p90_control_spread_scaled_treated               descriptive 1.84  
#>  smd_max_main_range                              descriptive 0.125 
#>  smd_max_main_sd                                 descriptive 0.0686
#>  smd_max_extended_range                          descriptive 0.215 
#>  smd_max_extended_sd                             descriptive 0.108 
#>  energy_distance_range                           descriptive 0.0398
#>  energy_distance_sd                              descriptive 0.0213
#>  ks_max_range                                    descriptive 0.222 
#>  ks_max_sd                                       descriptive 0.126 
#>  pair_recurrence_models_mean_treated             descriptive 0.4   
#>  pair_instability_models_mean_treated            exploratory 0.6   
#>  pair_recurrence_bootstrap_mean_treated          descriptive 0.62  
#>  treated_kept_share                              descriptive 1     
#>  treated_kept_share_min                          descriptive 1     
#>  bootstrap_rank_instability_treated              descriptive 0.0361
#>  reference_mean_absolute_rank_difference_treated descriptive 0.0766
#> 
#> Models (in_use = kept by the quality gate; for your own model, the row
#> `reference`, passed_quality_gate means within the balance limit):
#>               model_id                     algorithm logloss smd_max_main
#>  logistic_regression_1           logistic_regression   0.412        0.271
#>        random_forest_1                 random_forest   0.363        0.160
#>    gradient_boosting_1             gradient_boosting   0.351        0.285
#>              reference reference (glm, linear.logit)      NA        0.239
#>  smd_max_extended energy_distance ks_max treated_kept passed_quality_gate
#>             0.459          0.0953  0.357            1               FALSE
#>             0.244          0.0555  0.141            1                TRUE
#>             0.377          0.0885  0.135            1               FALSE
#>             0.239          0.0840  0.341            1                TRUE
#>  in_use
#>    TRUE
#>    TRUE
#>    TRUE
#>   FALSE
#> 
#> Gate notes:
#>   logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
#>       |SMD| main terms more than balance_margin above the best model's or
#>       above balance_cap; used anyway: fewer than two models passed
#>   random_forest_1: passed
#>   gradient_boosting_1: max |SMD| main terms more than balance_margin above the
#>       best model's or above balance_cap; used anyway: fewer than two models
#>       passed
#>   reference: your own model: an in-sample score from the matchit object, so
#>       its log-loss is not comparable and not computed; judged on balance only:
#>       within the balance limit; never one of the compared models
#> 
#> Quality gate: log-loss within 5 % of the best; relative balance rule: max
#>    |SMD| main terms at most 0.10 above the best model's and at most 0.25
#>    (limit here 0.25).
#> Matching: 1:1 nearest neighbour on the logit, with replacement, order
#>    "largest", caliper none; target: ATT.
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
# The default model set with one seed per tree algorithm (three models), to keep the example
# short; an analysis uses the default seeds = 1:3 (seven models, about a minute here).
propensity_sensitivity(matched, seeds = 1, bootstrap_samples = 5)
#> Warning: Fewer than two models passed the quality gate, so there is no set of similarly good models to compare. Every model is used below; read the values as a description of all models, not of a kept set.
#> Propensity model sensitivity
#> WARNING: fewer than two models passed the quality gate, so there is no set of
#> similarly good models to compare. Every model is used below; read the values
#> as a description of all models, not of a kept set.
#> 1. Quality gate: 3 of 3 models in use (3 algorithms: logistic regression,
#>    BART, GRF); log-loss within 5 % of the best; relative balance rule: max
#>    |SMD| main terms at most 0.10 above the best model's and at most 0.25
#>    (limit here 0.25).
#>    Fewer than two models passed, so every model is used. Failed:
#>      logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
#>        |SMD| main terms more than balance_margin above the best model's or
#>        above balance_cap
#>      bart_1: max |SMD| main terms more than balance_margin above the best
#>        model's or above balance_cap
#>      grf_probability_1: log-loss above (1 + quality_tolerance) x best; max
#>        |SMD| main terms more than balance_margin above the best model's or
#>        above balance_cap
#>    Treated patients kept by each model: logistic_regression_1 100 %, bart_1
#>      100 %, grf_probability_1 100 %
#> 2. Rank disagreement (primary; study mean over treated): 0.0669
#>    Equal weight per algorithm (exploratory): 0.0669
#>    Between/within-algorithm ratio (exploratory; median over treated): NA
#>    Pair instability (exploratory; study mean over treated): 0.611 of a
#>      possible 0.667 with 3 models. It is high whenever many untreated patients
#>      have similar scores, even when the models agree.
#>    Bootstrap rank instability of the logistic model (descriptive): 0.0361
#> 3. Balance spread (descriptive; range of max |SMD| main terms over the models
#>    in use): 0.145 (0.271 to 0.416)
#> Control spread (descriptive; mean scaled over treated): 1.22; 1 = the partners
#>    chosen by different models differ as much as a typical matched pair (scale:
#>    the matches of logistic_regression_1)
#> Your own model (glm, linear.logit): within the balance limit (judged on
#>    balance only; it is not one of the compared models); mean absolute rank
#>    difference from the models in use: 0.0747
#> Matching: 1:1 nearest neighbour on the logit, with replacement, order
#>    "largest", caliper none; target: ATT.
# }
```
