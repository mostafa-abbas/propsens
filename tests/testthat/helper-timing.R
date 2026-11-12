# Timing tests are opt-in. They run only when the environment variable
# PROPSENS_TIMING_TESTS is "true".
# Reason: elapsed time depends on the load from other work on the machine, so a slow run
# under load is noise, not a defect.

skip_unless_timing_tests <- function() {
  skip_if_not(identical(Sys.getenv("PROPSENS_TIMING_TESTS"), "true"),
              "timing tests are opt-in: set PROPSENS_TIMING_TESTS=true")
}

# Checks the elapsed seconds of a system.time() result against the limit. With
# PROPSENS_TIMING_LOG set to a log file it also writes one line per test to the log.
expect_elapsed_under <- function(timing, limit_seconds, label) {
  elapsed <- timing[["elapsed"]]
  log_file <- Sys.getenv("PROPSENS_TIMING_LOG")
  if (nzchar(log_file)) {
    cat(sprintf("timing test | %s | elapsed %.1f s | limit %g s | %s\n", label, elapsed,
                limit_seconds, if (elapsed < limit_seconds) "pass" else "FAIL"),
        file = log_file, append = TRUE)
  }
  expect_lt(elapsed, limit_seconds, label = paste(label, "elapsed seconds"))
}
