# Internal helpers shared by the four diagnostics.

# Models used for across-model statistics. With the gate, the models that passed the log-loss
# gate; if fewer than two passed, every model is used and the flag is TRUE. The oracle
# entry and the bootstrap block of a match set are never models. A non-NULL `models` (model ids)
# overrides the gate: the wrapper passes the final quality-gate set this way; the flag then still
# reports the outcome of the log-loss gate.
models_in_use <- function(x, use_gate = TRUE, models = NULL) {
  table <- x$models
  if (is.null(table) || !all(c("model_id", "algorithm", "passed_logloss_gate") %in%
                               names(table))) {
    stop("The object has no model table with `model_id`, `algorithm` and ",
         "`passed_logloss_gate`. Create it with fit_propensity_models() and ",
         "match_each_model().")
  }
  passed <- table$model_id[table$passed_logloss_gate]
  fewer_than_two <- length(passed) < 2
  if (!is.null(models)) {
    if (!is.character(models) || length(models) < 1 || anyNA(models) ||
          !all(models %in% table$model_id)) {
      stop("`models` must be a character vector of model ids of the object: ",
           paste(table$model_id, collapse = ", "), ".")
    }
    used <- table$model_id[table$model_id %in% models]
  } else {
    used <- if (!use_gate || fewer_than_two) table$model_id else passed
  }
  list(model_ids = used,
       algorithms = table$algorithm[match(used, table$model_id)],
       gate_applied = is.null(models) && use_gate && !fewer_than_two,
       models_supplied = !is.null(models),
       quality_gate_left_fewer_than_two_models = fewer_than_two)
}

check_true_false <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop("`", name, "` must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

check_match_set <- function(match_set) {
  if (!inherits(match_set, "match_set")) {
    stop("`match_set` must be a `match_set` from match_each_model().")
  }
  invisible(TRUE)
}

# Standard deviation of each row (n - 1 denominator); NA when there are fewer than two columns.
row_sd <- function(values) {
  if (ncol(values) < 2) return(rep(NA_real_, nrow(values)))
  sqrt(rowSums((values - rowMeans(values))^2) / (ncol(values) - 1))
}

# n x M integer matrix: row position of the control matched to each treated row, per model.
control_positions <- function(match_set, model_ids) {
  n <- length(match_set$treatment)
  positions <- matrix(NA_integer_, nrow = n, ncol = length(model_ids),
                      dimnames = list(NULL, model_ids))
  for (id in model_ids) {
    table <- match_set$matches[[id]]
    positions[match(table$treated_id, match_set$row_ids), id] <-
      match(table$control_id, match_set$row_ids)
  }
  positions
}

# Most frequent value; ties go to the smallest row position (the earliest row of the data).
modal_value <- function(values) {
  counts <- table(values)
  tied <- as.integer(names(counts)[counts == max(counts)])
  min(tied)
}

# Checks a covariate data frame and returns its numeric design matrix: factors and logicals
# become dummies (treatment contrasts), constant and linearly dependent columns are dropped
# with a message.
prepare_covariate_matrix <- function(covariates, n, caller) {
  if (!is.data.frame(covariates)) stop("`covariates` must be a data frame.")
  if (nrow(covariates) != n) {
    stop("`covariates` has ", nrow(covariates), " rows but the match set has ", n,
         " patients. Pass the covariates in the data's row order.")
  }
  character_columns <- names(covariates)[vapply(covariates, is.character, logical(1))]
  if (length(character_columns) > 0) {
    stop("Character column(s) in `covariates`: ", paste(character_columns, collapse = ", "),
         ". Convert them to factors first.")
  }
  if (anyNA(covariates)) stop("Missing values in `covariates`. Impute first.")
  constant <- names(covariates)[vapply(covariates, function(column) {
    length(unique(column)) < 2
  }, logical(1))]
  covariates <- covariates[, setdiff(names(covariates), constant), drop = FALSE]
  if (ncol(covariates) == 0) stop("Every column of `covariates` is constant.")
  design <- stats::model.matrix(~ ., data = covariates)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  constant_dummies <- colnames(design)[apply(design, 2, function(column) {
    length(unique(column)) < 2
  })]
  design <- design[, setdiff(colnames(design), constant_dummies), drop = FALSE]
  decomposition <- qr(scale(design))
  keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  dependent <- colnames(design)[-keep]
  design <- design[, keep, drop = FALSE]
  dropped <- c(constant, constant_dummies)
  if (length(dropped) > 0 || length(dependent) > 0) {
    message(caller, "(): dropped ",
            if (length(dropped) > 0) paste0("constant column(s): ",
                                            paste(dropped, collapse = ", ")),
            if (length(dropped) > 0 && length(dependent) > 0) "; ",
            if (length(dependent) > 0) paste0("linearly dependent column(s): ",
                                              paste(dependent, collapse = ", ")),
            ".")
  }
  list(design = design, dropped_constant = dropped, dropped_dependent = dependent)
}

# Rows transformed so that Euclidean distance equals Mahalanobis distance with covariance
# `cov(design)`: with cov = R'R (Cholesky), z = x R^-1.
whiten_rows <- function(design) {
  root <- chol(stats::cov(design))
  design %*% backsolve(root, diag(ncol(design)))
}
