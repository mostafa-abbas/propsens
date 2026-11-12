#' Simulate an example dataset for the propensity model sensitivity diagnostics
#'
#' Generates one dataset from simulated data-generating mechanisms (with the true propensity
#' and the outcome), for examples and for checking the diagnostics where the truth is known.
#' Six standard normal covariates `x1` to `x6`; in Scenario 5
#' four binary covariates `b1` to `b4` are added, `b4` rare at 5 %. The true treatment logit is
#' `alpha + X beta + nonlinear(X)`, where `alpha` is solved numerically so that the mean true
#' propensity over the generated patients equals `treated_share`. The outcome is
#' `y = mu0 + tau * treat + noise` with `mu0 = X gamma + lambda * nonlinear(X)`.
#'
#' Draw order: the covariates, the treatment and the outcome noise are drawn
#' first and do not depend on `lambda` or `heterogeneity`, so for a given scenario and seed
#' they are identical for every `lambda`; `lambda` enters only `mu0` and `y`.
#'
#' Scenarios:
#' * `"S1"` correct model: no nonlinear term.
#' * `"S2"` squared terms: `0.7 (x1^2 - 1) + 0.7 (x2^2 - 1)`.
#' * `"S3"` interactions: `0.9 x1 x2 + 0.9 x2 x3`.
#' * `"S4"` poor overlap: no nonlinear term, `beta` scaled by 2.5.
#' * `"S5"` mixed: S2 and S3 terms plus the four binary covariates (fixed coefficients).
#'
#' The function sets the random seed with `set.seed(seed)` and restores the previous state of
#' the random number generator on exit.
#'
#' @param n Number of patients (positive integer).
#' @param scenario One of `"S1"`, `"S2"`, `"S3"`, `"S4"`, `"S5"`.
#' @param lambda Strength of the nonlinear terms in the outcome model: `0` means the nonlinear
#'   terms affect treatment but not the outcome (harmless misspecification); larger values, such
#'   as `1.5`, make the misspecification harmful for the effect estimate.
#' @param treated_share Target mean true propensity, in (0, 1). Default `0.20`.
#' @param seed Integer seed.
#' @param heterogeneity If `TRUE`, the treatment effect is `tau = 1 + 0.8 * x1`; otherwise
#'   `tau = 1` for everyone.
#' @return A data frame with one row per patient and columns `id`, `x1` to `x6`, (`b1` to `b4`
#'   for Scenario 5), `treat` (0/1), `y`, `mu0`, `tau`, `true_logit`, `true_propensity`.
#'   Attributes `scenario`, `lambda`, `treated_share`, `seed`, `heterogeneity` and `alpha`
#'   record the settings.
#' @examples
#' example_data <- simulate_example_data(n = 500, scenario = "S2", lambda = 1.5, seed = 1)
#' head(example_data)
#' mean(example_data$treat)
#' attr(example_data, "alpha")
#' @importFrom stats rnorm rbinom plogis uniroot
#' @export
simulate_example_data <- function(n, scenario = c("S1", "S2", "S3", "S4", "S5"), lambda = 0,
                                  treated_share = 0.20, seed = 1, heterogeneity = FALSE) {
  scenario <- match.arg(scenario)
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 2 || n != round(n)) {
    stop("`n` must be a single integer of at least 2.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda)) {
    stop("`lambda` must be a single number.")
  }
  if (!is.numeric(treated_share) || length(treated_share) != 1 || is.na(treated_share) ||
        treated_share <= 0 || treated_share >= 1) {
    stop("`treated_share` must be a single number strictly between 0 and 1.")
  }
  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed)) {
    stop("`seed` must be a single integer.")
  }
  if (!is.logical(heterogeneity) || length(heterogeneity) != 1 || is.na(heterogeneity)) {
    stop("`heterogeneity` must be TRUE or FALSE.")
  }
  n <- as.integer(n)

  # Restore the caller's random number generator state on exit.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = globalenv())
    # nolint next: object_name_linter. ".Random.seed" is the base R name.
    on.exit(assign(".Random.seed", previous_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)

  # ---- fixed coefficients (including those of the binary covariates) ----
  beta <- c(0.5, -0.5, 0.4, 0.4, -0.3, 0.3)
  gamma <- c(1.0, -0.8, 1.0, 0.6, -0.6, 0.4)
  if (scenario == "S4") beta <- 2.5 * beta
  binary_prevalence <- c(0.5, 0.3, 0.2, 0.05)
  beta_binary <- c(0.4, -0.4, 0.3, 0.8)
  gamma_binary <- c(0.6, -0.5, 0.4, 1.0)

  # ---- draws that must not depend on lambda or heterogeneity ----
  x <- matrix(rnorm(n * 6), nrow = n, ncol = 6, dimnames = list(NULL, paste0("x", 1:6)))
  has_binary <- scenario == "S5"
  if (has_binary) {
    b <- vapply(binary_prevalence, function(prevalence) rbinom(n, 1, prevalence), integer(n))
    colnames(b) <- paste0("b", 1:4)
  }
  squared_terms <- 0.7 * (x[, 1]^2 - 1) + 0.7 * (x[, 2]^2 - 1)
  interaction_terms <- 0.9 * x[, 1] * x[, 2] + 0.9 * x[, 2] * x[, 3]
  nonlinear <- switch(scenario,
    S1 = rep(0, n),
    S2 = squared_terms,
    S3 = interaction_terms,
    S4 = rep(0, n),
    S5 = squared_terms + interaction_terms
  )
  linear_part <- drop(x %*% beta) + nonlinear
  if (has_binary) linear_part <- linear_part + drop(b %*% beta_binary)
  alpha <- uniroot(
    function(a) mean(plogis(a + linear_part)) - treated_share,
    interval = c(-50, 50), tol = 1e-10
  )$root
  true_logit <- alpha + linear_part
  true_propensity <- plogis(true_logit)
  treat <- rbinom(n, 1, true_propensity)
  noise <- rnorm(n)

  # ---- quantities that depend on lambda and heterogeneity ----
  mu0 <- drop(x %*% gamma) + lambda * nonlinear
  if (has_binary) mu0 <- mu0 + drop(b %*% gamma_binary)
  tau <- if (heterogeneity) 1 + 0.8 * x[, 1] else rep(1, n)
  y <- mu0 + tau * treat + noise

  out <- data.frame(id = seq_len(n), x)
  if (has_binary) out <- cbind(out, b)
  out$treat <- treat
  out$y <- y
  out$mu0 <- mu0
  out$tau <- tau
  out$true_logit <- true_logit
  out$true_propensity <- true_propensity
  attr(out, "scenario") <- scenario
  attr(out, "lambda") <- lambda
  attr(out, "treated_share") <- treated_share
  attr(out, "seed") <- seed
  attr(out, "heterogeneity") <- heterogeneity
  attr(out, "alpha") <- alpha
  out
}
