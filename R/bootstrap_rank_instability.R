#' Bootstrap rank instability of the logistic propensity model
#'
#' Descriptive. The sampling noise of the analyst's logistic model: logistic regression with the
#' formula as written (including any interactions or transformations) is refitted on
#' `bootstrap_samples` nonparametric bootstrap samples of all patients, the
#' original patients are ranked by each refit's logit (percentile rank
#' `rank(logit, ties = "average") / n`, on a 0 to 1 scale), and for each patient the standard
#' deviation of the rank over the refits is taken. The study-level value is the mean over treated
#' patients. A bootstrap sample with fewer than 10 treated or 10 untreated patients is skipped (the
#' rule of [match_each_model()]).
#'
#' It uses the covariates and the treatment only, no outcome. Report it next to rank
#' disagreement: it shows how much the logistic model's ranks move from sampling alone. It
#' describes the analyst's logistic model only, whatever distance the `matchit` object used.
#' The refits' scores are not clipped by default, unlike the cross-fitted scores (see
#' `clip_bootstrap`). These diagnostics have been studied in the author's simulations and
#' benchmark analyses; a methods paper is in preparation. No peer-reviewed evaluation exists yet.
#'
#' @param x One of: a `propensity_sensitivity` object or a `match_set` from
#'   [match_each_model()] (the stored bootstrap ranks are used, nothing is refitted); a `matchit`
#'   object (its formula is used); or a two-sided formula such as `treat ~ x1 + x2`.
#' @param data The data frame. Required with a formula; with a `matchit` object it is looked up
#'   from the object's call when `NULL`. Not used for stored ranks.
#' @param bootstrap_samples Number of bootstrap refits (default 20). Not used for
#'   stored ranks.
#' @param cores Number of cores (default 1). Not used for stored ranks.
#' @param seed Integer seed from which every other seed is derived; the same seed gives the
#'   bootstrap samples of [match_each_model()]. Not used for stored ranks.
#' @param clip_bootstrap `FALSE` (default) or `TRUE`: clip the refits' scores to
#'   `[0.001, 0.999]`, the default `clip` of [fit_propensity_models()]. Not used for stored
#'   ranks.
#' @return An object of class `"bootstrap_rank_instability"`: a list with `patients` (one row per
#'   patient: `id`, `treated`, `rank_instability`), `summary` (`mean_treated`, `samples_used`,
#'   `samples_requested`, `samples_skipped`) and `settings`. `mean_treated` is `NA` when fewer
#'   than two bootstrap samples are usable.
#' @examples
#' data("lalonde", package = "MatchIt")
#' instability <- bootstrap_rank_instability(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, bootstrap_samples = 5
#' )
#' instability
#' instability$summary$mean_treated
#' @export
bootstrap_rank_instability <- function(x, data = NULL, bootstrap_samples = 20, cores = 1,
                                       seed = 1, clip_bootstrap = FALSE) {
  if (inherits(x, "propensity_sensitivity")) x <- x$components$match_set
  if (inherits(x, "match_set")) {
    if (is.null(x$bootstrap)) {
      stop("This match set has no bootstrap block; use `bootstrap_samples` of at least 2 when ",
           "matching.")
    }
    return(instability_from_ranks(x$bootstrap$rank, x$treatment, x$row_ids, "stored ranks"))
  }
  caller <- parent.frame()
  if (inherits(x, "matchit")) {
    if (is.null(data)) {
      data <- tryCatch(eval(x$call$data, caller), error = function(e) NULL)
      if (!is.data.frame(data)) {
        stop("Could not find the data of the matchit object; pass the data with `data =`.")
      }
    }
    formula <- stats::as.formula(x$formula)
  } else if (inherits(x, "formula")) {
    if (is.null(data)) stop("With a formula, pass the data with `data =`.")
    formula <- x
  } else {
    stop("`x` must be a propensity_sensitivity object, a match_set, a matchit object or a ",
         "formula.")
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  check_single_integer(bootstrap_samples, "bootstrap_samples", minimum = 2)
  check_single_integer(cores, "cores", minimum = 1)
  check_single_integer(seed, "seed")
  check_true_false(clip_bootstrap, "clip_bootstrap")
  data <- treatment_as_0_1(formula, data)
  prepared <- prepare_model_data(formula, data)
  row_ids <- rownames(data)
  # Only the ranks are needed, so the bootstrap samples are not matched.
  options <- list(replace = TRUE, ratio = 1L, m_order = "largest", caliper = NULL)
  bootstrap <- bootstrap_logistic_matches(formula, data, prepared$treatment, row_ids, options,
                                          as.integer(bootstrap_samples), seed, cores,
                                          clip = if (clip_bootstrap) 0.001, match = FALSE)
  instability_from_ranks(bootstrap$rank, prepared$treatment, row_ids, "refitted")
}

# Skipped samples are columns of NA in the n x B rank matrix.
instability_from_ranks <- function(rank, treatment, row_ids, source) {
  usable <- rank[, !is.na(rank[1, ]), drop = FALSE]
  treated <- treatment == 1
  by_patient <- rep(NA_real_, nrow(rank))
  mean_treated <- NA_real_
  if (ncol(usable) >= 2) {
    by_patient <- unname(apply(usable, 1, stats::sd))
    mean_treated <- mean(by_patient[treated])
  }
  structure(list(
    patients = data.frame(id = row_ids, treated = as.integer(treatment),
                          rank_instability = by_patient, stringsAsFactors = FALSE),
    summary = list(mean_treated = mean_treated, samples_used = ncol(usable),
                   samples_requested = ncol(rank), samples_skipped = ncol(rank) - ncol(usable)),
    settings = list(source = source)
  ), class = "bootstrap_rank_instability")
}

#' @rdname bootstrap_rank_instability
#' @param ... Not used.
#' @export
print.bootstrap_rank_instability <- function(x, ...) {
  s <- x$summary
  cat("Bootstrap rank instability of the logistic model (mean over treated): ",
      format_number(s$mean_treated), "\n", sep = "")
  cat("SD of the percentile rank over ", s$samples_used, " bootstrap refits",
      if (s$samples_skipped > 0) paste0(" (", s$samples_skipped, " sample(s) skipped)"), "\n",
      sep = "")
  invisible(x)
}
