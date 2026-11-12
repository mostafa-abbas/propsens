# fit_propensity_models() defaults to logistic regression, BART and GRF.
# Tests that do not need BART or GRF pass the faster set explicitly.
faster_set <- c("logistic_regression", "random_forest", "gradient_boosting")

# Setup runs on lalonde often keep a single algorithm or fall back to every model; the two R
# warnings of propensity_sensitivity() for these cases are tested in their own tests and are
# silenced here. Any other warning still reaches the test.
quiet_gate <- function(expression) {
  withCallingHandlers(expression, warning = function(w) {
    text <- conditionMessage(w)
    if (startsWith(text, "Fewer than two models passed the quality gate") ||
          grepl("are random seeds of one algorithm", text, fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# The printed output as one string, with the line breaks and indentation of the wrapped lines
# collapsed to single spaces.
printed_text <- function(x) {
  gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))
}

# propensity_sensitivity() with the two gate warnings silenced (see quiet_gate()).
# The call is evaluated in the caller's frame, so that a matchit object's data is found there.
sensitivity_quietly <- function(...) {
  call <- sys.call()
  call[[1]] <- quote(propensity_sensitivity)
  quiet_gate(eval(call, parent.frame()))
}
