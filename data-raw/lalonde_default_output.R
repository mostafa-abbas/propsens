# data-raw/lalonde_default_output.R
#
# Pre-computes the printed output of propensity_sensitivity() with the default model set on the
# lalonde data (MatchIt), for the article "Getting started" (vignettes/getting-started.Rmd), so
# that the article builds in seconds. The default set takes about a minute per run.
#
# Writes (into vignettes/): lalonde-default-print.txt (print() of one run, seed = 1) and
# lalonde-default-seeds.txt (print() of seed_sensitivity() with seed = 1, 2, 3); the R warnings
# of the run are repeated at the top of the printed output. lalonde is public; only printed summaries
# are written.
#
# Usage, from the package root:  Rscript data-raw/lalonde_default_output.R

pkgload::load_all(".", quiet = TRUE)  # or library(propsens)
options(width = 80)
data("lalonde", package = "MatchIt")
matched <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                            data = lalonde, method = "nearest", distance = "glm",
                            link = "linear.logit", replace = TRUE)

sensitivity <- suppressWarnings(propensity_sensitivity(matched, seed = 1))
writeLines(utils::capture.output(print(sensitivity)), "vignettes/lalonde-default-print.txt")

by_seed <- seed_sensitivity(matched, seeds = 1:3)
writeLines(utils::capture.output(print(by_seed)), "vignettes/lalonde-default-seeds.txt")
