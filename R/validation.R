# Shared input checks, treatment coding and the parallel task runner.

# ---- treatment coding ---------------------------------------------------------------------

# A logical or two-level factor treatment as 0/1, with MatchIt's rule: TRUE is treated; for a
# factor, the level "0" (or the smaller of two numeric levels) is untreated when both levels are
# numbers, otherwise the first level is untreated. Other types are returned unchanged.
coerce_treatment <- function(treatment) {
  if (is.logical(treatment)) return(as.integer(treatment))
  if (!is.factor(treatment)) return(treatment)
  values <- levels(droplevels(treatment))
  if (length(values) != 2) {
    stop("A factor treatment must have exactly two levels (it has ", length(values), ").")
  }
  numeric_values <- suppressWarnings(as.numeric(values))
  untreated <- if (!anyNA(numeric_values)) {
    if (0 %in% numeric_values) values[numeric_values == 0] else values[which.min(numeric_values)]
  } else {
    values[1]
  }
  as.integer(as.character(treatment) != untreated)
}

# When the left-hand side of `formula` is a column of `data` holding a logical or factor
# treatment, the column is replaced by its 0/1 coding (a copy; row names are kept), so that every
# later step (glm, MatchIt, the bootstrap refits) sees the same coding.
treatment_as_0_1 <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3 || !is.data.frame(data)) {
    return(data)
  }
  if (!is.name(formula[[2]])) return(data)
  name <- as.character(formula[[2]])
  if (!name %in% names(data)) return(data)
  column <- data[[name]]
  if (is.logical(column) || is.factor(column)) data[[name]] <- coerce_treatment(column)
  data
}

# ---- argument checks ------------------------------------------------------------------------

check_matching_options <- function(replace, m_order, caliper) {
  if (!is.logical(replace) || length(replace) != 1 || is.na(replace)) {
    stop("`replace` must be TRUE or FALSE.")
  }
  allowed_orders <- c("largest", "smallest", "data", "random", "closest", "farthest")
  if (!is.character(m_order) || length(m_order) != 1 || !m_order %in% allowed_orders) {
    stop("`m_order` must be one of ", paste0('"', allowed_orders, '"', collapse = ", "), ".")
  }
  if (!is.null(caliper)) {
    if (!is.numeric(caliper) || length(caliper) != 1 || is.na(caliper) || caliper <= 0) {
      stop("`caliper` must be NULL or a single positive number (standard deviations of the ",
           "logit).")
    }
  }
  invisible(TRUE)
}

check_seeds <- function(seeds) {
  if (!is.numeric(seeds) || length(seeds) < 1 || anyNA(seeds) ||
        !isTRUE(all(seeds == seq_along(seeds)))) {
    stop("`seeds` must be 1:k, for example 1:3 for three seeds per tree algorithm: only the ",
         "number of seeds is used; set the randomness with `seed`.")
  }
  invisible(TRUE)
}

# The checks of the covariate data frame that balance_spread() and control_spread() would
# otherwise raise after every model is fitted.
check_covariate_frame <- function(covariates) {
  if (!is.data.frame(covariates)) stop("`covariates` must be a data frame.")
  character_columns <- names(covariates)[vapply(covariates, is.character, logical(1))]
  if (length(character_columns) > 0) {
    stop("Character column(s) in `covariates`: ", paste(character_columns, collapse = ", "),
         ". Convert them to factors first.")
  }
  if (anyNA(covariates)) stop("Missing values in `covariates`. Impute first.")
  invisible(TRUE)
}

check_proportion_or_null <- function(value, name) {
  if (is.null(value)) return(invisible(TRUE))
  if (!is.numeric(value) || length(value) != 1 || is.na(value) || value < 0 || value > 1) {
    stop("`", name, "` must be NULL or a single number between 0 and 1.")
  }
  invisible(TRUE)
}

# ---- parallel tasks -------------------------------------------------------------------------

propsens_state <- new.env(parent = emptyenv())

# Runs `task_function` over `tasks`, in parallel with parallel::mclapply() on Unix when
# cores > 1, serially otherwise (with one message per session on Windows). A task whose worker
# was stopped (for example by the operating system when memory runs out) comes back as NULL;
# that and an error in a task stop the run with the name of the task.
run_tasks <- function(tasks, task_function, cores, task_names = NULL, what = "task") {
  if (cores > 1 && .Platform$OS.type != "unix" && !isTRUE(propsens_state$windows_noted)) {
    propsens_state$windows_noted <- TRUE
    message("propsens: `cores` > 1 has no effect on Windows; the tasks run one after another.")
  }
  use_parallel <- cores > 1 && .Platform$OS.type == "unix"
  if (!use_parallel) return(lapply(tasks, task_function))
  results <- parallel::mclapply(tasks, task_function, mc.cores = cores, mc.preschedule = FALSE,
                                mc.set.seed = FALSE)
  if (is.null(task_names)) task_names <- paste(what, seq_along(tasks))
  for (i in seq_along(results)) {
    if (is.null(results[[i]])) {
      stop("The ", what, " \"", task_names[i], "\" returned no result: its worker process was ",
           "probably stopped (for example, out of memory). Rerun with fewer `cores`.")
    }
    if (inherits(results[[i]], "try-error")) {
      stop("The ", what, " \"", task_names[i], "\" failed: ",
           conditionMessage(attr(results[[i]], "condition")))
    }
  }
  results
}
