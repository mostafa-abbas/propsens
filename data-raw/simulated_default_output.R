# data-raw/simulated_default_output.R
#
# Pre-computes the printed output of propensity_sensitivity() with the default model set on one
# simulated dataset (simulate_example_data(): Scenario 1, 1,000 patients, seed 1), for the README
# and the article "Getting started", so that they build in seconds. The data are simulated: the
# example shows what the output looks like when the gate keeps models from several algorithms.
#
# Writes (into vignettes/): simulated-default-print.txt (print() of one run, seed = 1) and
# simulated-default-seeds.txt (print() of seed_sensitivity() with seed = 1, 2, 3).
#
# Usage, from the package root:  Rscript data-raw/simulated_default_output.R

pkgload::load_all(".", quiet = TRUE)  # or library(propsens)
options(width = 80)
simulated <- simulate_example_data(n = 1000, scenario = "S1", seed = 1)
formula <- treat ~ x1 + x2 + x3 + x4 + x5 + x6

sensitivity <- propensity_sensitivity(formula, data = simulated)
writeLines(utils::capture.output(print(sensitivity)), "vignettes/simulated-default-print.txt")

by_seed <- seed_sensitivity(formula, seeds = 1:3, data = simulated)
writeLines(utils::capture.output(print(by_seed)), "vignettes/simulated-default-seeds.txt")
