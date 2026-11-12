#' Fit the propensity model set with cross-fitting
#'
#' Fits every model of a propensity model set, the tree algorithms once per seed. The default is
#' the same as in [propensity_sensitivity()]: logistic regression with main effects, BART and the
#' GRF probability forest (three seeds each), with the same settings. A faster set of logistic
#' regression, random forests and gradient boosting is available with
#' `algorithms = c("logistic_regression", "random_forest", "gradient_boosting")` (packages
#' 'randomForest' and 'gbm'). These diagnostics have been studied in the author's simulations and
#' benchmark analyses; a methods paper is in preparation. No peer-reviewed evaluation exists yet.
#' Scores are
#' cross-fitted by default: the fold assignment
#' is stratified by treatment, derived from `seed`, and shared by every model; for each fold the
#' model is fitted on the other folds and predicts the held-out fold. Probabilities are clipped
#' to `[clip, 1 - clip]`, turned into logits, and ranked (percentile rank =
#' `rank(logit, ties = "average") / n`).
#'
#' @section Fixed settings:
#' Random forest: `ntree = 500`, `mtry = floor(sqrt(p))`, `nodesize = 10`.
#' Gradient boosting: `distribution = "bernoulli"`, `interaction.depth = 3`, `shrinkage = 0.05`,
#' `bag.fraction = 0.8`; the number of trees is frozen: when `n_trees` is
#' `NULL` it is chosen once by 5-fold cross-validation on the full data (`gbm` with
#' `cv.folds = 5`, capped at 3,000 trees, derived seed) and then used unchanged for every seed
#' and every fold; when `n_trees` is given it is used as is.
#'
#' @section Seeds and parallel runs:
#' Per-task seeds are derived in the main process from `seed` and the model's position in the
#' list of models (the order of `algorithms`, then the seed index), passed to each worker and
#' applied with [withr::with_seed()], which restores the random number generator state
#' afterwards. Results are identical for `cores = 1` and `cores > 1`. Because the position
#' enters the seed, a model's random stream depends on which algorithms come before it:
#' `random_forest_1` is a different fit with and without logistic regression in `algorithms`.
#' The fold split also comes from `seed`. On small data the scores, and therefore the models
#' kept by the quality gate and every summary of [propensity_sensitivity()], can change
#' noticeably with `seed`; rerun with two or three values of `seed` and report all of them
#' ([seed_sensitivity()] does this).
#' Parallel runs use [parallel::mclapply()] on Unix; on Windows the run is serial (with one
#' message per session). A worker that is stopped (for example, out of memory) stops the run
#' with the name of its model.
#'
#' @section The "as MatchIt does it" mode:
#' With `cross_fit = FALSE` the scores are in-sample predictions for logistic regression and
#' gradient boosting and out-of-bag predictions for the random forest. This mirrors what
#' `MatchIt::matchit()` returns for `distance = "glm"`, `"gbm"` and `"randomforest"`, so that
#' applied users can compare with their own object. In-sample tree scores are overfit and
#' shrink overlap, which is why cross-fitting is the default.
#'
#' @section Quality gate:
#' The models table reports the cross-fitted log-loss, Brier score and c-statistic (rank-based
#' AUC). `passed_logloss_gate` is `TRUE` when the log-loss is at most
#' `(1 + quality_tolerance)` times the best model's log-loss. The balance part of the
#' quality gate needs the matches and is applied later by [propensity_sensitivity()]: by
#' default a model is kept when its largest post-match absolute standardized mean difference on
#' main terms is at most 0.10 above the best model's and at most 0.25 (`balance_rule`).
#'
#' @param formula A two-sided formula such as `treat ~ x1 + x2`. The left-hand side is the
#'   treatment: numeric coded 0/1, logical (`TRUE` = treated) or a factor with two levels (coded
#'   as MatchIt does: the level `"0"`, or the smaller of two numeric levels, is untreated;
#'   otherwise the first level is untreated). Factor covariates are expanded to dummies for
#'   the tree algorithms.
#' @param data A data frame holding every variable in `formula`. Missing values are an error:
#'   impute first.
#' @param algorithms Which algorithms to fit. The default, `c("logistic_regression", "bart",
#'   "grf_probability")`, is also the default of
#'   [propensity_sensitivity()]: `"bart"` (`dbarts::bart2`, default priors, 75 trees, 500
#'   burn-in and 500 kept draws in each of 4 chains, one thread; the score is the
#'   posterior-mean probability) and `"grf_probability"` (`grf::probability_forest`, honest
#'   splitting, 2,000 trees, one thread; out-of-bag predictions when `cross_fit = FALSE`), both
#'   on the dummy-coded design matrix. The faster set, `"random_forest"` and
#'   `"gradient_boosting"` with logistic regression, needs the packages 'randomForest' and
#'   'gbm'; it has been studied less than the default set.
#'   Also exploratory, logistic regression variants on the same design matrix:
#'   `"logistic_splines"` (fitted once; natural cubic splines, `splines::ns(df = 3)`, on every
#'   column with at least 5 distinct values in the training rows whose two interior knots are
#'   distinct and inside its range, the other columns linear),
#'   `"lasso"` (`glmnet::cv.glmnet`, binomial, lambda.min from 5 internal folds that differ by
#'   seed) and `"lasso_interactions"` (the same on the main effects plus all two-way products of
#'   the design columns); the lasso variants need the package 'glmnet'.
#' @param seeds `1:k`: `k` tree models per tree algorithm (default `1:3`). Only the number of
#'   seeds is used, so any other vector is an error; the randomness is set by `seed`. Logistic
#'   regression is fitted once.
#' @param folds Number of cross-fitting folds (at least 2).
#' @param cross_fit `TRUE` for out-of-fold scores (default); `FALSE` for the "as MatchIt does
#'   it" mode described above.
#' @param n_trees `NULL` to choose the gradient boosting tree count once by cross-validation,
#'   or a positive integer to use as is.
#' @param quality_tolerance Relative log-loss tolerance of the quality gate (default 0.05).
#' @param clip Probabilities are clipped to `[clip, 1 - clip]` before the logit.
#' @param cores Number of cores for the model tasks (default 1).
#' @param seed Integer seed from which every other seed is derived.
#' @return An object of class `"propensity_model_set"`: a list with `logit`, `rank` and
#'   `probability` (n x M matrices whose column names are the model ids, for example
#'   `"logistic_regression_1"`, `"random_forest_2"`), `models` (data frame with `model_id`,
#'   `algorithm`, `seed_index`, `task_seed`, `n_trees`, `logloss`, `brier`, `c_statistic`,
#'   `passed_logloss_gate`), `folds` (integer fold assignment), `treatment`, `formula`,
#'   `logistic_log` (cross-fitting only: one row per logistic model and fold with the
#'   prediction `path` (`"formula"`, `"design_unseen_levels"` or
#'   `"design_fallback_predict_error"`), `glm_not_converged` and `glm_fitted_0_or_1`; `NULL`
#'   without logistic regression), `settings` (every argument plus the fixed tree settings and
#'   the chosen tree count) and `call`.
#' @examples
#' data("lalonde", package = "MatchIt")
#' model_set <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1, n_trees = 100,
#'   algorithms = c("logistic_regression", "random_forest", "gradient_boosting")
#' )
#' model_set
#' head(model_set$rank)
#' \donttest{
#' # The default algorithms (logistic regression, BART, GRF) with one seed per tree algorithm,
#' # to keep the example short; an analysis uses the default seeds = 1:3 (seven models).
#' default_algorithms <- fit_propensity_models(
#'   treat ~ age + educ + race + married + nodegree + re74 + re75,
#'   data = lalonde, seeds = 1
#' )
#' default_algorithms$models
#' }
#' @importFrom withr with_seed
#' @importFrom stats model.frame model.response model.matrix na.pass binomial predict
#' @export
fit_propensity_models <- function(formula, data,
                                  algorithms = c("logistic_regression", "bart",
                                                 "grf_probability"),
                                  seeds = 1:3, folds = 5, cross_fit = TRUE, n_trees = NULL,
                                  quality_tolerance = 0.05, clip = 0.001, cores = 1,
                                  seed = 1) {
  call <- match.call()
  algorithms <- match.arg(algorithms, propensity_algorithms, several.ok = TRUE)
  check_optional_algorithms(algorithms)
  check_single_integer(seed, "seed")
  check_single_integer(folds, "folds", minimum = 2)
  check_single_integer(cores, "cores", minimum = 1)
  check_seeds(seeds)
  if (!is.logical(cross_fit) || length(cross_fit) != 1 || is.na(cross_fit)) {
    stop("`cross_fit` must be TRUE or FALSE.")
  }
  if (!is.null(n_trees)) check_single_integer(n_trees, "n_trees", minimum = 1)
  if (!is.numeric(quality_tolerance) || length(quality_tolerance) != 1 ||
        is.na(quality_tolerance) || quality_tolerance < 0) {
    stop("`quality_tolerance` must be a single number of at least 0.")
  }
  if (!is.numeric(clip) || length(clip) != 1 || is.na(clip) || clip <= 0 || clip >= 0.5) {
    stop("`clip` must be a single number strictly between 0 and 0.5.")
  }
  data <- treatment_as_0_1(formula, data)
  prepared <- prepare_model_data(formula, data)
  n <- prepared$n

  tree_settings <- list(
    random_forest = list(ntree = 500L, mtry = max(1L, floor(sqrt(prepared$p))), nodesize = 10L),
    gradient_boosting = list(interaction.depth = 3L, shrinkage = 0.05, bag.fraction = 0.8,
                             n_trees_cap = 3000L, cv_folds = 5L),
    # dbarts::bart2 defaults (priors, trees, burn-in, draws, chains), one
    # thread so that a seed gives the same draws; grf::probability_forest defaults (honest).
    bart = list(n.trees = 75L, n.burn = 500L, n.samples = 500L, n.chains = 4L, n.threads = 1L),
    grf_probability = list(num.trees = 2000L, honesty = TRUE, num.threads = 1L),
    # Exploratory: logistic regression variants.
    logistic_splines = list(df = 3L, minimum_distinct_values = 5L),
    lasso = list(nfolds = 5L, s = "lambda.min"),
    lasso_interactions = list(nfolds = 5L, s = "lambda.min")
  )
  if ("lasso_interactions" %in% algorithms) {
    # Built once from the covariates only (never the treatment), so no fold sees another's labels.
    prepared$interaction_design <- two_way_interaction_design(prepared$design)
  }

  fold_id <- assign_stratified_folds(prepared$treatment, folds, derive_seed(seed, 0))

  n_trees_chosen <- n_trees
  n_trees_source <- "supplied"
  if ("gradient_boosting" %in% algorithms && is.null(n_trees)) {
    n_trees_chosen <- choose_n_trees(prepared$design, prepared$treatment, derive_seed(seed, 1),
                                     tree_settings$gradient_boosting)
    n_trees_source <- "cross-validation"
  }

  specs <- build_model_specs(algorithms, seeds, seed, n_trees_chosen)
  fit_task <- function(spec) {
    fit_one_model(spec, prepared, data, formula, fold_id, folds, cross_fit, tree_settings)
  }
  model_ids <- vapply(specs, `[[`, "", "model_id")
  results <- run_tasks(specs, fit_task, cores, task_names = model_ids, what = "model fit")

  # per model and fold, the logistic prediction path and glm's warnings.
  logistic_log <- do.call(rbind, lapply(seq_along(results), function(i) {
    log <- attr(results[[i]], "logistic_log")
    if (is.null(log)) NULL else data.frame(model_id = model_ids[i], log, stringsAsFactors = FALSE)
  }))
  warn_logistic_folds(logistic_log)
  probability <- do.call(cbind, results)
  probability <- pmin(pmax(probability, clip), 1 - clip)
  dimnames(probability) <- list(NULL, model_ids)
  logit <- log(probability / (1 - probability))
  rank <- apply(logit, 2, function(column) rank(column, ties.method = "average") / n)
  dimnames(rank) <- dimnames(logit)

  quality <- t(vapply(model_ids, function(id) model_quality(probability[, id], prepared$treatment),
                      numeric(3)))
  models <- data.frame(
    model_id = model_ids,
    algorithm = vapply(specs, `[[`, "", "algorithm"),
    seed_index = vapply(specs, `[[`, integer(1), "seed_index"),
    task_seed = vapply(specs, `[[`, integer(1), "task_seed"),
    n_trees = vapply(specs, `[[`, integer(1), "n_trees"),
    logloss = quality[, "logloss"],
    brier = quality[, "brier"],
    c_statistic = quality[, "c_statistic"],
    stringsAsFactors = FALSE
  )
  models$passed_logloss_gate <- models$logloss <= (1 + quality_tolerance) * min(models$logloss)
  rownames(models) <- NULL

  structure(list(
    logit = logit, rank = rank, probability = probability, models = models,
    folds = fold_id, treatment = prepared$treatment, formula = formula,
    logistic_log = logistic_log,
    settings = list(algorithms = algorithms, seeds = seeds, folds = folds, cross_fit = cross_fit,
                    n_trees = n_trees_chosen, n_trees_source = n_trees_source,
                    quality_tolerance = quality_tolerance, clip = clip, cores = cores,
                    seed = seed, fold_seed = derive_seed(seed, 0), tree_settings = tree_settings,
                    rng = rng_settings(), covariate_names = colnames(prepared$design)),
    call = call
  ), class = "propensity_model_set")
}

#' @rdname fit_propensity_models
#' @param x A `propensity_model_set` object.
#' @param digits Digits for the quality numbers.
#' @param ... Not used.
#' @export
print.propensity_model_set <- function(x, digits = 4, ...) {
  n <- length(x$treatment)
  cat("Propensity model set: ", ncol(x$logit), " model(s), ", n, " patients, ",
      sum(x$treatment), " treated\n", sep = "")
  scoring <- if (isTRUE(x$settings$cross_fit)) {
    paste0(x$settings$folds, "-fold cross-fitted scores")
  } else {
    "in-sample / out-of-bag scores (as MatchIt does it)"
  }
  cat("Scores: ", scoring, "; probabilities clipped to [", x$settings$clip, ", ",
      1 - x$settings$clip, "]\n", sep = "")
  if (!is.null(x$settings$n_trees)) {
    cat("Gradient boosting trees: ", x$settings$n_trees, " (", x$settings$n_trees_source, ")\n",
        sep = "")
  }
  cat("Log-loss gate: within ", 100 * x$settings$quality_tolerance,
      " % of the best model (balance part applied after matching)\n", sep = "")
  table <- x$models[, c("model_id", "logloss", "brier", "c_statistic", "passed_logloss_gate")]
  table$logloss <- round(table$logloss, digits)
  table$brier <- round(table$brier, digits)
  table$c_statistic <- round(table$c_statistic, digits)
  print(table, row.names = FALSE)
  invisible(x)
}

# ---- internal helpers ------------------------------------------------------------------

check_single_integer <- function(value, name, minimum = NULL) {
  ok <- is.numeric(value) && length(value) == 1 && !is.na(value) && value == round(value)
  if (ok && !is.null(minimum)) ok <- value >= minimum
  if (!ok) {
    stop("`", name, "` must be a single integer",
         if (!is.null(minimum)) paste0(" of at least ", minimum), ".")
  }
  invisible(TRUE)
}

prepare_model_data <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("`formula` must be a two-sided formula such as `treat ~ x1 + x2`.")
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  model_frame <- model.frame(formula, data = data, na.action = na.pass)
  if (anyNA(model_frame)) {
    stop("Missing values in the treatment or the covariates. Impute first; ",
         "fit_propensity_models() does not handle missing values.")
  }
  treatment <- coerce_treatment(model.response(model_frame))
  if (!is.numeric(treatment) || !all(treatment %in% c(0, 1))) {
    stop("The treatment (left-hand side of `formula`) must be numeric coded 0/1 ",
         "(0 = untreated, 1 = treated), logical, or a factor with two levels. Recode it first.")
  }
  if (length(unique(treatment)) < 2) stop("The treatment must contain both 0 and 1.")
  design <- model.matrix(formula, model_frame)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (ncol(design) == 0) stop("`formula` must contain at least one covariate.")
  # Factor and character columns of the model frame (after factor() and other calls inside the
  # formula), for the check of levels absent from a fold or a bootstrap sample.
  covariate_frame <- model_frame[-1]
  level_columns <- covariate_frame[vapply(covariate_frame, function(column) {
    is.factor(column) || is.character(column)
  }, logical(1))]
  list(treatment = as.integer(treatment), design = design, n = nrow(design), p = ncol(design),
       level_columns = level_columns)
}

derive_seed <- function(seed, index) {
  as.integer((as.numeric(seed) * 1000003 + as.numeric(index) * 7919) %% 2147483647)
}

rng_settings <- function() {
  list(kind = "Mersenne-Twister", normal_kind = "Inversion", sample_kind = "Rejection")
}

with_task_seed <- function(task_seed, code) {
  rng <- rng_settings()
  withr::with_seed(task_seed, code, .rng_kind = rng$kind, .rng_normal_kind = rng$normal_kind,
                   .rng_sample_kind = rng$sample_kind)
}

assign_stratified_folds <- function(treatment, folds, fold_seed) {
  fold_id <- integer(length(treatment))
  with_task_seed(fold_seed, {
    for (group in c(0L, 1L)) {
      rows <- which(treatment == group)
      fold_id[rows] <- sample(rep(seq_len(folds), length.out = length(rows)))
    }
  })
  fold_id
}

build_model_specs <- function(algorithms, seeds, seed, n_trees) {
  specs <- list()
  index <- 0L
  for (algorithm in algorithms) {
    seed_indices <- if (algorithm %in% single_fit_algorithms) 1L else seq_along(seeds)
    for (seed_index in seed_indices) {
      index <- index + 1L
      specs[[index]] <- list(
        model_id = paste0(algorithm, "_", seed_index),
        algorithm = algorithm,
        seed_index = as.integer(seed_index),
        task_seed = derive_seed(seed, 100L + index),
        n_trees = if (algorithm == "gradient_boosting") as.integer(n_trees) else NA_integer_
      )
    }
  }
  specs
}

choose_n_trees <- function(design, treatment, n_trees_seed, gbm_settings) {
  frame <- data.frame(treatment = treatment, design, check.names = TRUE)
  # gbm's cross-validation attaches gbm, which prints a start-up message; it is suppressed.
  best <- with_task_seed(n_trees_seed, suppressPackageStartupMessages({
    invisible(utils::capture.output(
      fit <- gbm::gbm(treatment ~ ., data = frame, distribution = "bernoulli",
                      n.trees = gbm_settings$n_trees_cap,
                      interaction.depth = gbm_settings$interaction.depth,
                      shrinkage = gbm_settings$shrinkage, bag.fraction = gbm_settings$bag.fraction,
                      cv.folds = gbm_settings$cv_folds, n.cores = 1, verbose = FALSE)
    ))
    gbm::gbm.perf(fit, method = "cv", plot.it = FALSE)
  }))
  best <- as.integer(best)
  if (best >= gbm_settings$n_trees_cap) {
    warning("The cross-validated tree count hit the cap of ", gbm_settings$n_trees_cap,
            " trees; consider supplying `n_trees` after your own tuning.")
  }
  best
}

fit_one_model <- function(spec, prepared, data, formula, fold_id, folds, cross_fit,
                          tree_settings) {
  n <- prepared$n
  with_task_seed(spec$task_seed, {
    if (cross_fit) {
      probability <- numeric(n)
      fold_log <- list()
      for (k in seq_len(folds)) {
        test <- fold_id == k
        if (!any(test)) next
        held_out <- predict_held_out(spec, prepared, data, formula, !test, test, tree_settings)
        if (!is.null(attr(held_out, "logistic_log"))) {
          fold_log[[length(fold_log) + 1]] <- data.frame(fold = k, attr(held_out, "logistic_log"))
        }
        probability[test] <- held_out
      }
      if (length(fold_log) > 0) attr(probability, "logistic_log") <- do.call(rbind, fold_log)
      probability
    } else {
      predict_as_matchit(spec, prepared, data, formula, tree_settings)
    }
  })
}

predict_held_out <- function(spec, prepared, data, formula, train, test, tree_settings) {
  switch(spec$algorithm,
    logistic_regression = {
      # each fold records which path predicted it and whether glm warned of
      # non-convergence or of fitted probabilities of 0 or 1 (warnings still pass through).
      warned <- character(0)
      log_row <- function(path) {
        data.frame(path = path,
                   glm_not_converged = any(grepl("did not converge", warned, fixed = TRUE)),
                   glm_fitted_0_or_1 = any(grepl("fitted probabilities numerically 0 or 1",
                                                 warned, fixed = TRUE)),
                   stringsAsFactors = FALSE)
      }
      note <- function(w) {
        warned <<- c(warned, conditionMessage(w))
        if (is_known_glm_warning(w)) invokeRestart("muffleWarning")
      }
      if (has_unseen_levels(prepared, train, test)) {
        probability <- withCallingHandlers(logistic_on_design(prepared, train, test),
                                           warning = note)
        attr(probability, "logistic_log") <- log_row("design_unseen_levels")
        return(probability)
      }
      train_data <- data[train, , drop = FALSE]
      fit <- withCallingHandlers(stats::glm(formula, data = train_data, family = binomial()),
                                 warning = note)
      # a rank-deficient fit that did not converge (many dummy columns, separation) can
      # leave non-finite values in its QR, and predict() then stops; the design-matrix path
      # gives the same linear predictor with non-estimable coefficients set to 0.
      path <- "formula"
      fallback <- function(e) {
        path <<- "design_fallback_predict_error"
        withCallingHandlers(logistic_on_design(prepared, train, test), warning = note)
      }
      test_data <- data[test, , drop = FALSE]
      probability <- tryCatch(as.numeric(predict(fit, newdata = test_data, type = "response")),
                              error = fallback)
      attr(probability, "logistic_log") <- log_row(path)
      probability
    },
    random_forest = {
      settings <- tree_settings$random_forest
      fit <- randomForest::randomForest(
        x = prepared$design[train, , drop = FALSE],
        y = factor(prepared$treatment[train], levels = c(0, 1)),
        ntree = settings$ntree, mtry = settings$mtry, nodesize = settings$nodesize
      )
      as.numeric(predict(fit, newdata = prepared$design[test, , drop = FALSE],
                         type = "prob")[, "1"])
    },
    gradient_boosting = {
      fit <- fit_gbm(prepared$design[train, , drop = FALSE], prepared$treatment[train],
                     spec$n_trees, tree_settings$gradient_boosting)
      as.numeric(predict(fit, newdata = as.data.frame(prepared$design[test, , drop = FALSE]),
                         n.trees = spec$n_trees, type = "response"))
    },
    bart = fit_bart(prepared$design[train, , drop = FALSE], prepared$treatment[train],
                    prepared$design[test, , drop = FALSE], tree_settings$bart),
    grf_probability = fit_grf_probability(prepared$design[train, , drop = FALSE],
                                          prepared$treatment[train],
                                          prepared$design[test, , drop = FALSE],
                                          tree_settings$grf_probability),
    logistic_splines = fit_logistic_splines(prepared$design[train, , drop = FALSE],
                                            prepared$treatment[train],
                                            prepared$design[test, , drop = FALSE],
                                            tree_settings$logistic_splines),
    lasso = fit_lasso(prepared$design[train, , drop = FALSE], prepared$treatment[train],
                      prepared$design[test, , drop = FALSE], tree_settings$lasso),
    lasso_interactions = fit_lasso(prepared$interaction_design[train, , drop = FALSE],
                                   prepared$treatment[train],
                                   prepared$interaction_design[test, , drop = FALSE],
                                   tree_settings$lasso_interactions)
  )
}

# TRUE when a factor or character column of the model frame takes, in the held-out rows, a level
# that the training rows lack (rare levels and cross-fitting). glm() drops such levels and
# predict() then stops; logistic_on_design() is used instead. Otherwise the formula path is kept
# unchanged. The columns come from the model frame, so a factor made inside the formula
# (`factor(x)`) is checked as well.
has_unseen_levels <- function(prepared, train, test) {
  any(vapply(prepared$level_columns, function(column) {
    !all(unique(column[test]) %in% unique(column[train]))
  }, logical(1)))
}

# The two glm() warnings that the logistic log records per fold or bootstrap sample.
is_known_glm_warning <- function(w) {
  text <- conditionMessage(w)
  grepl("did not converge", text, fixed = TRUE) ||
    grepl("fitted probabilities numerically 0 or 1", text, fixed = TRUE)
}

# One warning for all folds in which glm() warned, instead of one per fold.
warn_logistic_folds <- function(logistic_log) {
  if (is.null(logistic_log)) return(invisible(FALSE))
  warned <- logistic_log$glm_not_converged | logistic_log$glm_fitted_0_or_1
  if (!any(warned)) return(invisible(FALSE))
  kinds <- c(if (any(logistic_log$glm_not_converged)) "did not converge",
             if (any(logistic_log$glm_fitted_0_or_1)) "fitted probabilities of 0 or 1")
  warning("Logistic regression: glm() warned in ", sum(warned), " of ", length(warned),
          " cross-fitting folds (", paste(kinds, collapse = "; "), "). This usually means that ",
          "some covariate patterns predict treatment almost perfectly (separation). The folds ",
          "are listed in `logistic_log` of the model set.", call. = FALSE)
  invisible(TRUE)
}

# Logistic regression on the full-data design matrix (the dummy columns of every level, as for
# the tree models). A level absent from the training rows gives an all-zero column whose
# coefficient is not estimable; it is set to 0 for prediction, so held-out patients with that
# level are scored as the reference level of that factor.
logistic_on_design <- function(prepared, train, test) {
  x_train <- cbind(1, prepared$design[train, , drop = FALSE])
  fit <- stats::glm.fit(x_train, prepared$treatment[train], family = stats::binomial())
  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  as.numeric(stats::plogis(cbind(1, prepared$design[test, , drop = FALSE]) %*% coefficients))
}

predict_as_matchit <- function(spec, prepared, data, formula, tree_settings) {
  switch(spec$algorithm,
    logistic_regression = {
      fit <- stats::glm(formula, data = data, family = binomial())
      as.numeric(predict(fit, type = "response"))
    },
    random_forest = {
      settings <- tree_settings$random_forest
      fit <- randomForest::randomForest(
        x = prepared$design, y = factor(prepared$treatment, levels = c(0, 1)),
        ntree = settings$ntree, mtry = settings$mtry, nodesize = settings$nodesize
      )
      as.numeric(predict(fit, type = "prob")[, "1"])   # out-of-bag, as MatchIt
    },
    gradient_boosting = {
      fit <- fit_gbm(prepared$design, prepared$treatment, spec$n_trees,
                     tree_settings$gradient_boosting)
      as.numeric(predict(fit, newdata = as.data.frame(prepared$design), n.trees = spec$n_trees,
                         type = "response"))
    },
    bart = fit_bart(prepared$design, prepared$treatment, NULL, tree_settings$bart),
    grf_probability = fit_grf_probability(prepared$design, prepared$treatment, NULL,
                                          tree_settings$grf_probability),
    logistic_splines = fit_logistic_splines(prepared$design, prepared$treatment, NULL,
                                            tree_settings$logistic_splines),
    lasso = fit_lasso(prepared$design, prepared$treatment, NULL, tree_settings$lasso),
    lasso_interactions = fit_lasso(prepared$interaction_design, prepared$treatment, NULL,
                                   tree_settings$lasso_interactions)
  )
}

fit_gbm <- function(design, treatment, n_trees, gbm_settings) {
  gbm::gbm.fit(x = as.data.frame(design), y = treatment, distribution = "bernoulli",
               n.trees = n_trees, interaction.depth = gbm_settings$interaction.depth,
               shrinkage = gbm_settings$shrinkage, bag.fraction = gbm_settings$bag.fraction,
               keep.data = FALSE, verbose = FALSE)
}

compute_c_statistic <- function(score, treatment) {
  n1 <- sum(treatment == 1)
  n0 <- sum(treatment == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[treatment == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

model_quality <- function(probability, treatment) {
  c(logloss = -mean(treatment * log(probability) + (1 - treatment) * log(1 - probability)),
    brier = mean((probability - treatment)^2),
    c_statistic = compute_c_statistic(probability, treatment))
}

# ---- BART and the GRF probability forest --------------------------------------
# Each fit draws its own integer seed from the task's random-number
# stream, so a given task seed gives the same scores on one core or many. Causal forests are not
# used: they need the outcome, and the design is outcome-free.
propensity_algorithms <- c("logistic_regression", "random_forest", "gradient_boosting", "bart",
                           "grf_probability", "logistic_splines", "lasso", "lasso_interactions")
# Fitted once (no seed enters the fit); the others once per seed.
single_fit_algorithms <- c("logistic_regression", "logistic_splines")

check_optional_algorithms <- function(algorithms) {
  packages <- c(bart = "dbarts", grf_probability = "grf", lasso = "glmnet",
                lasso_interactions = "glmnet", random_forest = "randomForest",
                gradient_boosting = "gbm")
  needed <- packages[intersect(algorithms, names(packages))]
  for (algorithm in names(needed)) {
    # gbm prints a start-up message ("no longer under development") when it is loaded.
    available <- suppressPackageStartupMessages(requireNamespace(needed[[algorithm]],
                                                                 quietly = TRUE))
    if (!available) {
      stop("Algorithm \"", algorithm, "\" needs the package '", needed[[algorithm]], "'.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

draw_library_seed <- function() sample.int(.Machine$integer.max, 1L)

# Posterior-mean probability of treatment from dbarts::bart2 (probit for a 0/1 outcome). With
# x_test = NULL the in-sample posterior mean is returned.
fit_bart <- function(x_train, treatment, x_test, settings) {
  library_seed <- draw_library_seed()
  fit <- dbarts::bart2(x_train, treatment, test = x_test, n.trees = settings$n.trees,
                       n.burn = settings$n.burn, n.samples = settings$n.samples,
                       n.chains = settings$n.chains, n.threads = settings$n.threads,
                       seed = library_seed, verbose = FALSE)
  as.numeric(stats::fitted(fit, type = "ev", sample = if (is.null(x_test)) "train" else "test"))
}

# grf::probability_forest on the dummy-coded design (honest splitting by default). With
# x_test = NULL the out-of-bag predictions are returned, as for the random forest.
fit_grf_probability <- function(x_train, treatment, x_test, settings) {
  library_seed <- draw_library_seed()
  fit <- grf::probability_forest(x_train, factor(treatment, levels = c(0, 1)),
                                 num.trees = settings$num.trees, honesty = settings$honesty,
                                 num.threads = settings$num.threads, seed = library_seed)
  predictions <- if (is.null(x_test)) predict(fit)$predictions else
    predict(fit, newdata = x_test)$predictions
  as.numeric(predictions[, "1"])
}

# ---- Exploratory: logistic regression variants ---------------------------------------------------
# Logistic regression with natural cubic splines (splines::ns, `df` degrees of freedom) on every
# design column with at least `minimum_distinct_values` distinct values in the training rows whose
# interior knots (ns()'s quantiles) are distinct and strictly inside the range; other columns
# (dummies, binary, few-valued or heavily concentrated covariates) enter linearly. The spline
# basis (its knots) comes from the training rows and is applied to the held-out rows. Coefficients
# that are not estimable are set to 0, as in logistic_on_design().
fit_logistic_splines <- function(x_train, treatment, x_test, settings) {
  splined <- vapply(seq_len(ncol(x_train)), function(j) {
    spline_eligible(x_train[, j], settings$df, settings$minimum_distinct_values)
  }, logical(1))
  expand <- function(x, bases) {
    parts <- lapply(seq_len(ncol(x)), function(j) {
      if (splined[j]) stats::predict(bases[[j]], x[, j]) else x[, j, drop = FALSE]
    })
    do.call(cbind, parts)
  }
  bases <- lapply(seq_len(ncol(x_train)), function(j) {
    if (splined[j]) splines::ns(x_train[, j], df = settings$df) else NULL
  })
  train_matrix <- cbind(1, expand(x_train, bases))
  fit <- stats::glm.fit(train_matrix, treatment, family = stats::binomial())
  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  new_matrix <- if (is.null(x_test)) train_matrix else cbind(1, expand(x_test, bases))
  as.numeric(stats::plogis(new_matrix %*% coefficients))
}

spline_eligible <- function(values, df, minimum_distinct_values) {
  if (length(unique(values)) < minimum_distinct_values) return(FALSE)
  knots <- stats::quantile(values, seq_len(df - 1) / df, names = FALSE)
  all(diff(c(min(values), knots, max(values))) > 0)
}

# glmnet lasso (binomial), lambda chosen by internal cross-validation on the training rows; the
# internal folds are drawn from the task's random-number stream, so each seed gives other folds.
fit_lasso <- function(x_train, treatment, x_test, settings) {
  internal_folds <- sample(rep_len(seq_len(settings$nfolds), nrow(x_train)))
  fit <- glmnet::cv.glmnet(x_train, treatment, family = "binomial", foldid = internal_folds,
                           alpha = 1)
  new_x <- if (is.null(x_test)) x_train else x_test
  as.numeric(stats::predict(fit, newx = new_x, s = settings$s, type = "response"))
}

# Main effects plus all two-way products of the design columns (no squares).
two_way_interaction_design <- function(design) {
  pairs <- utils::combn(ncol(design), 2)
  products <- design[, pairs[1, ], drop = FALSE] * design[, pairs[2, ], drop = FALSE]
  colnames(products) <- paste(colnames(design)[pairs[1, ]], colnames(design)[pairs[2, ]],
                              sep = ":")
  cbind(design, products)
}
