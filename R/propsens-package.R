#' propsens: Propensity Model Sensitivity Diagnostics for Matched Studies
#'
#' Outcome-free diagnostics that measure how much a matched study depends on
#' the choice of propensity model. The definitions are given in the help pages and in
#' `vignette("definitions", package = "propsens")`.
#'
#' The `importFrom` lines below declare the packages the diagnostics will use,
#' so that `R CMD check` sees every entry of `Imports` used.
#'
#' @importFrom MatchIt matchit
#' @importFrom cobalt bal.tab bal.init bal.compute
#' @importFrom randomForest randomForest
#' @importFrom gbm gbm
#' @importFrom stats glm
#' @importFrom parallel mclapply
#' @keywords internal
"_PACKAGE"
