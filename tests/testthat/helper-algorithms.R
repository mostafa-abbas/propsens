# fit_propensity_models() defaults to the paper's set (logistic regression, BART, GRF).
# Tests that do not need BART or GRF pass the faster set explicitly.
faster_set <- c("logistic_regression", "random_forest", "gradient_boosting")
