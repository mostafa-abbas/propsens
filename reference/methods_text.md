# Methods and Results sentences for a propensity model sensitivity result

Builds the sentences to report a
[`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
result, from the result itself: the model set, the quality gate and the
models it kept (by algorithm), the matching settings, the share of
treated patients kept, the seed, and the main summaries. After a
fallback (fewer than two models passed the gate) the text says that all
models were used and describes them as all models. Edit the sentences to
fit your paper; call the gate criteria pre-specified only if you fixed
them before running the check.

## Usage

``` r
methods_text(x, digits = 3)
```

## Arguments

- x:

  A `propensity_sensitivity` object.

- digits:

  Decimal places for the rank summaries.

## Value

A named character vector with `methods` and `results`, one paragraph
each.

## Examples

``` r
data("lalonde", package = "MatchIt")
sensitivity <- propensity_sensitivity(
  treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
  seeds = 1, n_trees = 50, bootstrap_samples = 5, replace = TRUE
)
#> Warning: Fewer than two models passed the quality gate, so there is no set of similarly good models to compare. Every model is used below; read the values as a description of all models, not of a kept set.
cat(methods_text(sensitivity), sep = "\n\n")
#> We described the sensitivity of the matched design to the propensity model with the R package propsens (version 0.1.0). Propensity models from logistic regression, random forest and gradient boosting (1 seed for each tree method; 3 models; 5-fold cross-fitted scores; seed 1) were each used for 1:1 nearest-neighbour matching with replacement. Models were kept if their cross-fitted log-loss was within 5 % of the best model's and their largest post-match standardized mean difference on main terms was within 0.10 of the best model's and at most 0.25; fewer than two of the 3 models met these criteria, so all models were used, and all treated patients were matched. For all models we report the mean rank disagreement over treated patients, the bootstrap rank instability of the logistic model, and the range of the largest standardized mean difference on main terms. No outcome data were used. These summaries are descriptive; they have not been validated as predictors of bias.
#> 
#> Mean rank disagreement over treated patients was 0.069 across all models; the bootstrap rank instability of the logistic model was 0.036. As exploratory summaries, pair instability was 0.63 (of a possible 0.67 with 3 models) and the median between/within-algorithm ratio was not available. The largest standardized mean difference on main terms ranged from 0.16 to 0.54 across all models.
```
