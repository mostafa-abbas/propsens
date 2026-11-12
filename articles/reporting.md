# Reading the plots and reporting in a paper

The check belongs in the design section of a matched study, next to the
balance table, and it is done before any outcome is looked at. This
vignette shows how to read the plots and how to build the Methods and
Results text from the result, so that every number and setting comes
from the run.

``` r

library(propsens)
data("lalonde", package = "MatchIt")
matched <- MatchIt::matchit(
  treat ~ age + educ + race + married + nodegree + re74 + re75,
  data = lalonde, distance = "glm", link = "linear.logit", replace = TRUE
)
# The faster model set and small settings so the vignette builds quickly; use the defaults
# (logistic regression, BART and the GRF probability forest) in a real analysis.
sensitivity <- propensity_sensitivity(
  matched, algorithms = c("logistic_regression", "random_forest", "gradient_boosting"),
  seeds = 1:2, n_trees = 100, bootstrap_samples = 5
)
```

## Reading the three plots

**Patients plot** (`plot(x, type = "patients")`). Each point is a
treated patient: right = the kept models place this patient differently
among all patients; up = the kept models give this patient different
kinds of untreated partners. The plot describes how the disagreement is
spread over the treated patients; it is a description, and individual
points should not be read as flags for individual patients. Points near
the bottom left mean that these models rank the patient similarly and
choose similar partners. Colour shows whether the patient’s logistic
score is outside the range of the untreated patients’ scores.

**Models plot** (`plot(x, type = "models")`). Each row is a model, each
panel a balance statistic, each axis starting at 0. If the points of a
panel line up, the balance you report does not depend on the model; if
they spread out, it does: say so, and show this plot. Filled points are
kept by the quality gate, hollow points are not; your own model has its
own colour and is judged on balance only. The extended-terms panel is
usually driven by sparse squares and interactions; report the main-terms
panel.

**Recurrence plot** (`plot(x, type = "recurrence")`). The per-patient
pair recurrence behind pair instability (1 minus its mean over treated
patients across the kept models; exploratory). Values near the floor (1
/ number of kept models) are expected when untreated patients are dense,
even when the models agree; report pair instability next to rank
disagreement, never alone.

## What to write in Methods and Results

[`methods_text()`](https://mostafa-abbas.github.io/propsens/reference/methods_text.md)
builds both paragraphs from the result: the model set, the gate and how
many models it kept (by algorithm), the matching settings, the share of
treated patients kept, the seed and the main summaries. After a fallback
it says that all models were used. Edit the sentences to fit your paper.

``` r

text <- methods_text(sensitivity)
```

**Methods:**

> We described the sensitivity of the matched design to the propensity
> model with the R package propsens (version 0.1.0). Propensity models
> from logistic regression, random forest and gradient boosting (2 seeds
> for each tree method; 5 models; 5-fold cross-fitted scores; seed 1)
> were each used for 1:1 nearest-neighbour matching with replacement.
> Models were kept if their cross-fitted log-loss was within 5 % of the
> best model’s and their largest post-match standardized mean difference
> on main terms was within 0.10 of the best model’s and at most 0.25; 3
> of 5 models were kept (random forest: 2; gradient boosting: 1), and
> all treated patients were matched. For the kept models we report the
> mean rank disagreement over treated patients, the bootstrap rank
> instability of the logistic model, and the range of the largest
> standardized mean difference on main terms. No outcome data were used.
> These summaries are descriptive; they have not been validated as
> predictors of bias.

**Results:**

> Mean rank disagreement over treated patients was 0.029 across the kept
> models; the bootstrap rank instability of the logistic model was
> 0.036. As exploratory summaries, pair instability was 0.52 (of a
> possible 0.67 with 3 models) and the median between/within-algorithm
> ratio was 4.5. The largest standardized mean difference on main terms
> ranged from 0.16 to 0.20 across the kept models.

If you reran the check with several values of `seed`
([`seed_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/seed_sensitivity.md)),
add: “The check was repeated with seeds 1, 2 and 3; all results are
reported.”

The text states the balance rule in effect and, after a fallback, that
all models were used. Call the criteria pre-specified only if you fixed
them before running the check.

## If you judge the disagreement to be large

No thresholds exist yet, so this is a judgement you state and justify.

1.  Show the models plot (balance by model).
2.  Choose the design by balance, never by the outcome.
3.  Decide before looking at outcomes whether you will report the effect
    estimate under each kept model, and then report all of them.
4.  Look at the colour in the patients plot: is the disagreement
    concentrated among treated patients whose logistic score lies
    outside the range of the untreated patients’ scores (outside the
    support, an overlap problem)? Describe this for the study as a
    whole.

## What the check does not do

It does not detect unmeasured confounding, it does not fix a design, and
it does not tell you which model is right. When logistic regression is
correct, the tree models can still disagree with it: the check describes
that the *choice* matters, not which choice is correct. It cannot tell
whether a disagreement matters for the effect estimate, and a small
value is not evidence that the estimate is unbiased. It covers 1:1
nearest-neighbour matching for the ATT only.
