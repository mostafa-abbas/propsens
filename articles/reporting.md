# Reading the plots and reporting in a paper

The check belongs in the design section of a matched study, next to the
balance table, and it is done before any outcome is looked at. This
vignette shows how to read the output and what to write.

``` r

library(propsens)
data("lalonde", package = "MatchIt")
matched <- MatchIt::matchit(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, distance = "glm", link = "linear.logit", replace = TRUE
)
# An earlier, faster model set and small settings so the vignette builds quickly; use the
# defaults (logistic regression, BART and the GRF probability forest) in a real analysis.
sensitivity <- propensity_sensitivity(
  matched, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
  seeds = 1:2, n_trees = 100, bootstrap_samples = 5
)
numbers <- sensitivity$summary$table
value <- function(name) numbers$with_gate[numbers$quantity == name]
```

## Reading the three plots

**Patients plot** (`plot(x, type = "patients")`). Each point is a
treated patient. Right = the models disagree about where this patient
sits among all patients. Up = the models give this patient different
kinds of untreated partners. The upper right is what makes a design
fragile. A cloud that hugs the bottom left means the design is robust to
the choice of model. Colour shows whether the patient’s logistic score
is outside the range of the untreated patients’ scores: fragile patients
outside the support point to an overlap problem; fragile patients inside
the support point to a region with few comparable untreated patients.

**Models plot** (`plot(x, type = "models")`). Each row is a model, each
panel a balance statistic. If the points of a panel line up, the balance
you report does not depend on the model. If they spread out, it does:
say so, and show this plot. Hollow points fail the quality gate; your
own model has its own colour and is judged on balance only.

**Recurrence plot** (`plot(x, type = "recurrence")`). The per-patient
pair recurrence behind pair instability (1 minus its mean over treated
patients across models), the second, exploratory summary printed under
rank disagreement. Values near the floor (1 / number of models) are
common when untreated patients are dense; report pair instability next
to rank disagreement, not alone.

## What to write in Methods

> We assessed sensitivity of the matched design to the propensity model
> with the `propsens` package: logistic regression, Bayesian additive
> regression trees and a generalized random forest (fixed settings,
> 5-fold cross-fitted scores; seven models in all) were each used for
> 1:1 nearest-neighbour matching with replacement; we report the mean
> rank disagreement with the bootstrap rank instability of the logistic
> model, the balance spread across models and, as descriptions, pair
> instability and the control spread, all computed without outcomes.

## What to write in Results

The numbers below come from the run above, not typed by hand:

> Mean rank disagreement was 0.029 (study level); the bootstrap rank
> instability of the logistic model was 0.036. Exploratory summaries:
> pair instability 0.524 and a between/within ratio of 4.5. As a
> description only, 59 % of treated patients received untreated partners
> that differed by more than one typical match distance across models.
> The maximum standardized mean difference (including squares and
> interactions) ranged from 0.23 to 0.24 across models (Figure S1).

If the quality gate left fewer than two models, add: “Fewer than two
propensity models met the pre-specified quality criteria (log-loss
within 5 % of the best model’s; maximum standardized mean difference on
main terms at most 0.10 above the best model’s and at most 0.25), so all
models were used.” With `balance_rule = "absolute"`, describe that rule
instead; [`print()`](https://rdrr.io/r/base/print.html) states the rule
in effect.

## When the result says “fragile”

1.  Show the models plot (balance by model).
2.  Choose the design by balance, never by the outcome.
3.  Report the range of the estimates across the models in a supplement.
4.  Look at the treated patients with high control spread: do they share
    a region of the covariates with few untreated patients? Say so in
    the Discussion.

## What the check does not do

It does not detect unmeasured confounding, it does not fix a fragile
design, and it does not tell you which model is right. When logistic
regression is correct, the tree models can still disagree with it: the
check flags that the *choice* matters, not which choice is correct. It
covers 1:1 nearest-neighbour matching for the ATT only.
