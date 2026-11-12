# Fit the propensity model set with cross-fitting

Fits every model of a propensity model set, the tree algorithms once per
seed. The default is the same as in
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md):
logistic regression with main effects, BART and the GRF probability
forest (three seeds each), with the same settings. A faster set of
logistic regression, random forests and gradient boosting is available
with
`algorithms = c("logistic_regression", "random_forest", "gradient_boosting")`
(packages 'randomForest' and 'gbm'). These diagnostics have been studied
in the author's simulations and benchmark analyses; a methods paper is
in preparation. No peer-reviewed evaluation exists yet. Scores are
cross-fitted by default: the fold assignment is stratified by treatment,
derived from `seed`, and shared by every model; for each fold the model
is fitted on the other folds and predicts the held-out fold.
Probabilities are clipped to `[clip, 1 - clip]`, turned into logits, and
ranked (percentile rank = `rank(logit, ties = "average") / n`).

## Usage

``` r
fit_propensity_models(
  formula,
  data,
  algorithms = c("logistic_regression", "bart", "grf_probability"),
  seeds = 1:3,
  folds = 5,
  cross_fit = TRUE,
  n_trees = NULL,
  quality_tolerance = 0.05,
  clip = 0.001,
  cores = 1,
  seed = 1
)

# S3 method for class 'propensity_model_set'
print(x, digits = 4, ...)
```

## Arguments

- formula:

  A two-sided formula such as `treat ~ x1 + x2`. The left-hand side is
  the treatment: numeric coded 0/1, logical (`TRUE` = treated) or a
  factor with two levels (coded as MatchIt does: the level `"0"`, or the
  smaller of two numeric levels, is untreated; otherwise the first level
  is untreated). Factor covariates are expanded to dummies for the tree
  algorithms.

- data:

  A data frame holding every variable in `formula`. Missing values are
  an error: impute first.

- algorithms:

  Which algorithms to fit. The default,
  `c("logistic_regression", "bart", "grf_probability")`, is also the
  default of
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md):
  `"bart"` ([`dbarts::bart2`](https://rdrr.io/pkg/dbarts/man/bart.html),
  default priors, 75 trees, 500 burn-in and 500 kept draws in each of 4
  chains, one thread; the score is the posterior-mean probability) and
  `"grf_probability"`
  ([`grf::probability_forest`](https://rdrr.io/pkg/grf/man/probability_forest.html),
  honest splitting, 2,000 trees, one thread; out-of-bag predictions when
  `cross_fit = FALSE`), both on the dummy-coded design matrix. The
  faster set, `"random_forest"` and `"gradient_boosting"` with logistic
  regression, needs the packages 'randomForest' and 'gbm'; it has been
  studied less than the default set. Also exploratory, logistic
  regression variants on the same design matrix: `"logistic_splines"`
  (fitted once; natural cubic splines, `splines::ns(df = 3)`, on every
  column with at least 5 distinct values in the training rows whose two
  interior knots are distinct and inside its range, the other columns
  linear), `"lasso"`
  ([`glmnet::cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
  binomial, lambda.min from 5 internal folds that differ by seed) and
  `"lasso_interactions"` (the same on the main effects plus all two-way
  products of the design columns); the lasso variants need the package
  'glmnet'.

- seeds:

  `1:k`: `k` tree models per tree algorithm (default `1:3`). Only the
  number of seeds is used, so any other vector is an error; the
  randomness is set by `seed`. Logistic regression is fitted once.

- folds:

  Number of cross-fitting folds (at least 2).

- cross_fit:

  `TRUE` for out-of-fold scores (default); `FALSE` for the "as MatchIt
  does it" mode described above.

- n_trees:

  `NULL` to choose the gradient boosting tree count once by
  cross-validation, or a positive integer to use as is.

- quality_tolerance:

  Relative log-loss tolerance of the quality gate (default 0.05).

- clip:

  Probabilities are clipped to `[clip, 1 - clip]` before the logit.

- cores:

  Number of cores for the model tasks (default 1).

- seed:

  Integer seed from which every other seed is derived.

- x:

  A `propensity_model_set` object.

- digits:

  Digits for the quality numbers.

- ...:

  Not used.

## Value

An object of class `"propensity_model_set"`: a list with `logit`, `rank`
and `probability` (n x M matrices whose column names are the model ids,
for example `"logistic_regression_1"`, `"random_forest_2"`), `models`
(data frame with `model_id`, `algorithm`, `seed_index`, `task_seed`,
`n_trees`, `logloss`, `brier`, `c_statistic`, `passed_logloss_gate`),
`folds` (integer fold assignment), `treatment`, `formula`,
`logistic_log` (cross-fitting only: one row per logistic model and fold
with the prediction `path` (`"formula"`, `"design_unseen_levels"` or
`"design_fallback_predict_error"`), `glm_not_converged` and
`glm_fitted_0_or_1`; `NULL` without logistic regression), `settings`
(every argument plus the fixed tree settings and the chosen tree count)
and `call`.

## Fixed settings

Random forest: `ntree = 500`, `mtry = floor(sqrt(p))`, `nodesize = 10`.
Gradient boosting: `distribution = "bernoulli"`,
`interaction.depth = 3`, `shrinkage = 0.05`, `bag.fraction = 0.8`; the
number of trees is frozen: when `n_trees` is `NULL` it is chosen once by
5-fold cross-validation on the full data (`gbm` with `cv.folds = 5`,
capped at 3,000 trees, derived seed) and then used unchanged for every
seed and every fold; when `n_trees` is given it is used as is.

## Seeds and parallel runs

Per-task seeds are derived in the main process from `seed` and the
model's position in the list of models (the order of `algorithms`, then
the seed index), passed to each worker and applied with
[`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
which restores the random number generator state afterwards. Results are
identical for `cores = 1` and `cores > 1`. Because the position enters
the seed, a model's random stream depends on which algorithms come
before it: `random_forest_1` is a different fit with and without
logistic regression in `algorithms`. The fold split also comes from
`seed`. On small data the scores, and therefore the models kept by the
quality gate and every summary of
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md),
can change noticeably with `seed`; rerun with two or three values of
`seed` and report all of them
([`seed_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/seed_sensitivity.md)
does this). Parallel runs use
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
Unix; on Windows the run is serial (with one message per session). A
worker that is stopped (for example, out of memory) stops the run with
the name of its model.

## The "as MatchIt does it" mode

With `cross_fit = FALSE` the scores are in-sample predictions for
logistic regression and gradient boosting and out-of-bag predictions for
the random forest. This mirrors what
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
returns for `distance = "glm"`, `"gbm"` and `"randomforest"`, so that
applied users can compare with their own object. In-sample tree scores
are overfit and shrink overlap, which is why cross-fitting is the
default.

## Quality gate

The models table reports the cross-fitted log-loss, Brier score and
c-statistic (rank-based AUC). `passed_logloss_gate` is `TRUE` when the
log-loss is at most `(1 + quality_tolerance)` times the best model's
log-loss. The balance part of the quality gate needs the matches and is
applied later by
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md):
by default a model is kept when its largest post-match absolute
standardized mean difference on main terms is at most 0.10 above the
best model's and at most 0.25 (`balance_rule`).

## Examples

``` r
data("lalonde", package = "MatchIt")
model_set <- fit_propensity_models(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, seeds = 1, n_trees = 100,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
)
model_set
#> Propensity model set: 3 model(s), 614 patients, 185 treated
#> Scores: 5-fold cross-fitted scores; probabilities clipped to [0.001, 0.999]
#> Gradient boosting trees: 100 (supplied)
#> Log-loss gate: within 5 % of the best model (balance part applied after matching)
#>               model_id logloss  brier c_statistic passed_logloss_gate
#>  logistic_regression_1  0.4123 0.1338      0.8612               FALSE
#>        random_forest_1  0.3634 0.1143      0.9014                TRUE
#>    gradient_boosting_1  0.3513 0.1110      0.9091                TRUE
head(model_set$rank)
#>      logistic_regression_1 random_forest_1 gradient_boosting_1
#> [1,]             0.7117264       0.7719870           0.8061889
#> [2,]             0.5781759       0.7622150           0.8192182
#> [3,]             0.7508143       0.9568404           0.9299674
#> [4,]             0.9706840       0.9918567           0.9983713
#> [5,]             0.8583062       0.9291531           0.9714984
#> [6,]             0.8631922       0.9649837           0.9820847
# \donttest{
# The default algorithms (logistic regression, BART, GRF) with one seed per tree algorithm,
# to keep the example short; an analysis uses the default seeds = 1:3 (seven models).
default_algorithms <- fit_propensity_models(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, seeds = 1
)
default_algorithms$models
#>                model_id           algorithm seed_index task_seed n_trees
#> 1 logistic_regression_1 logistic_regression          1   1799822      NA
#> 2                bart_1                bart          1   1807741      NA
#> 3     grf_probability_1     grf_probability          1   1815660      NA
#>     logloss     brier c_statistic passed_logloss_gate
#> 1 0.4122854 0.1337668   0.8612424               FALSE
#> 2 0.3529691 0.1112218   0.9048762                TRUE
#> 3 0.3777964 0.1208646   0.8924463               FALSE
# }
```
