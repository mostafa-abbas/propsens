# Spell check of the help pages and vignettes against inst/WORDLIST (Language: en-GB). It runs
# only when the 'spelling' package is installed, and is skipped on CRAN.
if (requireNamespace("spelling", quietly = TRUE)) {
  spelling::spell_check_test(vignettes = TRUE, error = FALSE, skip_on_cran = TRUE)
}
