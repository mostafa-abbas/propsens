# Hand-built objects for the diagnostic tests. They carry only the fields the diagnostics read.

fake_model_set <- function(logit, algorithms, treatment, rank = NULL, passed = TRUE,
                           formula = NULL) {
  n <- nrow(logit)
  seed_index <- stats::ave(seq_along(algorithms), algorithms, FUN = seq_along)
  ids <- paste0(algorithms, "_", seed_index)
  colnames(logit) <- ids
  if (is.null(rank)) rank <- apply(logit, 2, function(column) rank(column) / n)
  rank <- matrix(rank, nrow = n, dimnames = list(NULL, ids))
  structure(list(
    logit = logit, rank = rank, probability = stats::plogis(logit),
    models = data.frame(model_id = ids, algorithm = algorithms,
                        seed_index = as.integer(seed_index),
                        passed_logloss_gate = rep(passed, length.out = length(ids)),
                        stringsAsFactors = FALSE),
    treatment = as.integer(treatment), formula = formula
  ), class = "propensity_model_set")
}

# `pairs`: named list, one entry per model, each a two-column matrix (treated row, control row).
fake_match_set <- function(treatment, pairs, algorithms = "random_forest", bootstrap = NULL) {
  n <- length(treatment)
  row_ids <- as.character(seq_len(n))
  ids <- names(pairs)
  matches <- lapply(pairs, function(p) {
    data.frame(treated_id = row_ids[p[, 1]], control_id = row_ids[p[, 2]], weight = 1,
               pair_id = seq_len(nrow(p)), stringsAsFactors = FALSE)
  })
  weights <- vapply(pairs, function(p) {
    w <- numeric(n)
    w[p[, 1]] <- 1
    counts <- table(p[, 2])
    w[as.integer(names(counts))] <- as.numeric(counts)
    w
  }, numeric(n))
  dimnames(weights) <- list(row_ids, ids)
  structure(list(
    matches = matches, weights = weights,
    treated_kept = vapply(pairs, function(p) nrow(p) / sum(treatment), numeric(1)),
    oracle = NULL, bootstrap = bootstrap, model_ids = ids,
    models = data.frame(model_id = ids, algorithm = rep(algorithms, length.out = length(ids)),
                        seed_index = seq_along(ids), passed_logloss_gate = TRUE,
                        stringsAsFactors = FALSE),
    treatment = as.integer(treatment), row_ids = row_ids
  ), class = "match_set")
}
