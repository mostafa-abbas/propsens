# Changelog

## propsens 0.1.0

First public release.

- [`propensity_sensitivity()`](https://mostafa-abbas.github.io/propsens/reference/propensity_sensitivity.md)
  checks, before any outcome is used, how much a 1:1 matched design
  depends on the choice of propensity model. It takes a formula and data
  with your matching settings, or an existing `MatchIt` object, whose
  settings it reads.
- The model set is logistic regression, BART and the GRF probability
  forest (three seeds for each tree algorithm: seven models), with
  cross-fitted scores; every model is matched with the same settings.
- A quality gate keeps only the models that predict treatment about as
  well as the best model and whose matched groups are about as well
  balanced, and warns when fewer than two pass.
- For the kept models it reports rank disagreement (how differently they
  order the treated patients) and pair instability (how often they
  choose a different untreated partner), with the bootstrap rank
  instability of the logistic model and the spread of balance across
  models.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
  results, and the articles explain the method step by step with the
  public right heart catheterization study.
