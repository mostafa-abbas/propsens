# Simulate an example dataset for the propensity model sensitivity diagnostics

Generates one dataset from the data-generating mechanisms of the
simulation study. Six standard normal covariates `x1` to `x6`; in
Scenario 5 four binary covariates `b1` to `b4` are added, `b4` rare at 5
%. The true treatment logit is `alpha + X beta + nonlinear(X)`, where
`alpha` is solved numerically so that the mean true propensity over the
generated patients equals `treated_share`. The outcome is
`y = mu0 + tau * treat + noise` with
`mu0 = X gamma + lambda * nonlinear(X)`.

## Usage

``` r
simulate_example_data(
  n,
  scenario = c("S1", "S2", "S3", "S4", "S5"),
  lambda = 0,
  treated_share = 0.2,
  seed = 1,
  heterogeneity = FALSE
)
```

## Arguments

- n:

  Number of patients (positive integer).

- scenario:

  One of `"S1"`, `"S2"`, `"S3"`, `"S4"`, `"S5"`.

- lambda:

  Strength of the nonlinear terms in the outcome model: `0` means the
  nonlinear terms affect treatment but not the outcome (harmless
  misspecification); the plan uses `1.5` for the harmful case.

- treated_share:

  Target mean true propensity, in (0, 1). Default `0.20`.

- seed:

  Integer seed.

- heterogeneity:

  If `TRUE`, the treatment effect is `tau = 1 + 0.8 * x1`; otherwise
  `tau = 1` for everyone.

## Value

A data frame with one row per patient and columns `id`, `x1` to `x6`,
(`b1` to `b4` for Scenario 5), `treat` (0/1), `y`, `mu0`, `tau`,
`true_logit`, `true_propensity`. Attributes `scenario`, `lambda`,
`treated_share`, `seed`, `heterogeneity` and `alpha` record the
settings.

## Details

Draw order: the covariates, the treatment and the outcome noise are
drawn first and do not depend on `lambda` or `heterogeneity`, so for a
given scenario and seed they are identical for every `lambda`; `lambda`
enters only `mu0` and `y`.

Scenarios:

- `"S1"` correct model: no nonlinear term.

- `"S2"` squared terms: `0.7 (x1^2 - 1) + 0.7 (x2^2 - 1)`.

- `"S3"` interactions: `0.9 x1 x2 + 0.9 x2 x3`.

- `"S4"` poor overlap: no nonlinear term, `beta` scaled by 2.5.

- `"S5"` mixed: S2 and S3 terms plus the four binary covariates (fixed
  coefficients).

The function sets the random seed with `set.seed(seed)` and restores the
previous state of the random number generator on exit.

## Examples

``` r
example_data <- simulate_example_data(n = 500, scenario = "S2", lambda = 1.5, seed = 1)
head(example_data)
#>   id         x1          x2          x3         x4          x5         x6 treat
#> 1  1 -0.6264538  0.07730312  1.13496509  0.8500435 -0.88614959 -1.8054836     0
#> 2  2  0.1836433 -0.29686864  1.11193185 -0.9253130 -1.92225490 -0.6780407     0
#> 3  3 -0.8356286 -1.18324224 -0.87077763  0.8935812  1.61970074 -0.4733581     0
#> 4  4  1.5952808  0.01129269  0.21073159 -0.9410097  0.51926990  1.0274171     0
#> 5  5  0.3295078  0.99160104  0.06939565  0.5389521 -0.05584993 -0.5973876     0
#> 6  6 -0.8204684  1.59396745 -1.66264885 -0.1819744  0.69641761  1.1598494     0
#>            y        mu0 tau true_logit true_propensity
#> 1 -0.2739991 -0.9154677   1  -3.017625      0.04663598
#> 2  1.5431866 -0.1120324   1  -2.689377      0.06360314
#> 3 -2.1412086 -1.2815731   1  -2.438997      0.08024694
#> 4  2.5214621  1.9040781   1  -1.028983      0.26328140
#> 5  0.1058257 -1.2300108   1  -2.948719      0.04979711
#> 6 -3.0807998 -2.5467910   1  -3.019151      0.04656817
mean(example_data$treat)
#> [1] 0.218
attr(example_data, "alpha")
#> [1] -2.062843
```
