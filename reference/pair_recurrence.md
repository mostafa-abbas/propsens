# Pair recurrence across models and across bootstrap resamples

A descriptive index of how often the same untreated partner is chosen,
reported next to pair instability (1 minus the recurrence across
models). The across-models and bootstrap versions are kept separate and
are never mixed.

## Usage

``` r
pair_recurrence(match_set, use_gate = TRUE, models = NULL)

# S3 method for class 'pair_recurrence'
print(x, ...)
```

## Arguments

- match_set:

  A `match_set` from
  [`match_each_model()`](https://mostafa-abbas.github.io/propsens/reference/match_each_model.md).

- use_gate:

  `TRUE` (default) to use the models that passed the log-loss gate;
  `FALSE` to use every model. The bootstrap version does not depend on
  the gate.

- models:

  Optional character vector of model ids to use instead of a gate; it
  overrides `use_gate`. The wrapper
  [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  passes the final quality-gate set (log-loss and balance) this way.

- x:

  A `pair_recurrence` object.

- ...:

  Not used.

## Value

An object of class `"pair_recurrence"`: a list with `patients` (one row
per patient: `row`, `id`, `treated`, `modal_control_models`,
`pair_recurrence_models`, `modal_control_bootstrap`,
`pair_recurrence_bootstrap`, `pair_recurrence_bootstrap_uncorrected`;
`NA` for untreated patients), `summary` (means over treated,
`floor_models`, `ceiling_models`, `ceiling_bootstrap`,
`expected_ceiling_bootstrap_uncorrected`, `bootstrap_samples_used`),
`models_used`, `quality_gate_left_fewer_than_two_models` and `settings`.

## Details

**Across models.** For treated patient i, the modal partner is the
untreated partner chosen most often by the models in use; the recurrence
is the share of those models choosing it (see the denominator rule
below). Floor `1 / M` (every model picks a different untreated partner),
ceiling 1.

**Across bootstrap resamples** (when the match set has a bootstrap
block: logistic regression only, so it measures sampling instability).
The modal partner is taken over the resamples. The corrected recurrence
divides the number of resamples in which i is matched to its modal
partner by the number of resamples in which **both** i and that partner
are present; its ceiling is 1. The uncorrected version divides by the
number of resamples in which i is present; its expected ceiling is about
`1 - 1/e = 0.632`, because a given partner is missing from about 37 % of
bootstrap samples. Skipped resamples are ignored.

Expect the across-models version near its floor whenever untreated
patients are dense: the index measures crowding more than model
dependence.

## Tie rule and denominator rule

- **Ties.** When two or more untreated partners are chosen equally
  often, the modal partner is the one that comes first in the data
  (smallest row position). The same rule is used across models and
  across bootstrap resamples. With three models that all choose
  different untreated partners, the recurrence is 1/3 whichever partner
  is called modal.

- **Denominator across models.** The share is taken over the models in
  use that matched the patient. Without a caliper every model matches
  every treated patient, so the denominator is M, the number of models
  in use. With a caliper a model may leave a patient unmatched; that
  model is left out of the denominator, and a patient unmatched by every
  model gets `NA`. The reported floor is `1 / M`.

- **Denominator across bootstrap resamples.** Corrected: resamples in
  which both the patient and its modal partner are present. Uncorrected:
  resamples in which the patient is present. Skipped resamples are left
  out of both.

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
  data = lalonde, seeds = 1, n_trees = 100,
  algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
)
match_set <- match_each_model(model_set, lalonde, replace = TRUE, bootstrap_samples = 20)
recurrence <- pair_recurrence(match_set)
recurrence
#> Pair recurrence: 2 model(s) in use (after the log-loss gate)
#> Across models: mean over treated 0.5811 (floor 0.5, ceiling 1)
#> Across 20 bootstrap resamples: corrected 0.3119 (ceiling 1); uncorrected 0.2067 (expected ceiling 0.6321)
```
