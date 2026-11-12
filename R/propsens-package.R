#' propsens: Propensity Model Sensitivity Diagnostics for Matched Studies
#'
#' Outcome-free diagnostics of how much a matched design depends on which of several similarly
#' good propensity models is used. The definitions are given in the help pages and in
#' `vignette("definitions", package = "propsens")`. These diagnostics have been studied in the
#' author's simulations and benchmark analyses; a methods paper is in preparation. No
#' peer-reviewed evaluation exists yet.
#'
#' The `importFrom` lines below declare packages of `Imports` that are otherwise used only with
#' `::`, so that `R CMD check` sees every entry used.
#'
#' @importFrom MatchIt matchit
#' @importFrom cobalt bal.tab bal.init bal.compute
#' @importFrom stats glm
#' @importFrom parallel mclapply
#' @importFrom utils combn
#' @keywords internal
"_PACKAGE"
