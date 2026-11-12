# data-raw/heart_catheterization_example.R
#
# Pre-computes the figures and the tables of the heart catheterization example used in the
# articles "propsens at a glance" (vignettes/at-a-glance.Rmd) and "How propsens works, step by
# step" (vignettes/how-it-works.Rmd), so that the articles build in seconds.
#
# Data, all public:
# - Covariates and treatment: ATbounds::RHC, the right heart catheterization study (Connors et
#   al., JAMA 1996; 5,735 patients, 2,184 catheterized within 24 hours). Its dummy columns are
#   turned back into the 50 covariates of the original file (factors for the categorical ones);
#   the two derived columns surv2md1 and wt0 are not used.
# - Outcome: death within 30 days (dth30). ATbounds::RHC does not contain it (its `survival`
#   column is survival to the end of follow-up), so it is read from the public raw file of the
#   Vanderbilt Department of Biostatistics, https://hbiostat.org/data/repo/rhc.csv, after its
#   SHA-256 is checked. Only dth30 is taken from that file.
#
# Writes (into vignettes/): four PNG figures (rhc-*.png), rhc-models.csv (one row per model:
# fit, balance, treated kept, gate, matched odds ratio) and rhc-study.csv (study-level numbers). Aggregates only; no patient rows,
# scores or weights are written.
#
# The propensity models take about 16 minutes on one core. The fitted result is cached in
# data-raw/cache/ (not part of the package or of the repository); delete the cache to refit.
#
# Usage, from the package root:  Rscript data-raw/heart_catheterization_example.R

pkgload::load_all(".", quiet = TRUE)  # or library(propsens)
library(ggplot2)

rhc_url <- "https://hbiostat.org/data/repo/rhc.csv"
rhc_sha256 <- "9ef4ab578be4b40ad5d97d3a7e08ffdc1f9f76aeeefee51b4996e4221556f8e8"
cache_path <- "data-raw/cache/rhc_sensitivity_seed1.rds"
figure_dir <- "vignettes"

# ---- covariates and treatment from ATbounds ------------------------------------------------
rebuild_rhc_covariates <- function() {
  rhc <- get(utils::data("RHC", package = "ATbounds", envir = environment()))
  from_dummies <- function(levels_by_column, reference) {
    out <- rep(reference, nrow(rhc))
    for (column in names(levels_by_column)) out[rhc[[column]] == 1] <- levels_by_column[[column]]
    out
  }
  yes_no <- function(column) ifelse(rhc[[column]] == 1, "Yes", "No")
  cohort <- data.frame(
    swang1 = rhc$RHC,
    cat1 = from_dummies(c(cat1_CHF = "CHF", cat1_Cirrhosis = "Cirrhosis",
                          cat1_Colon_Cancer = "Colon Cancer", cat1_Coma = "Coma",
                          cat1_COPD = "COPD", cat1_Lung_Cancer = "Lung Cancer",
                          cat1_MOSF_Malignancy = "MOSF w/Malignancy",
                          cat1_MOSF_Sepsis = "MOSF w/Sepsis"), "ARF"),
    cat2 = from_dummies(c(cat2_Cirrhosis = "Cirrhosis", cat2_Colon_Cancer = "Colon Cancer",
                          cat2_Coma = "Coma", cat2_Lung_Cancer = "Lung Cancer",
                          cat2_MOSF_Malignancy = "MOSF w/Malignancy",
                          cat2_MOSF_Sepsis = "MOSF w/Sepsis"), "None"),
    ca = from_dummies(c(ca_Metastatic = "Metastatic", ca_Yes = "Yes"), "No"),
    stringsAsFactors = FALSE)
  for (v in c("cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx", "liverhx",
              "gibledhx", "malighx", "immunhx", "transhx", "amihx", "age")) cohort[[v]] <- rhc[[v]]
  cohort$sex <- ifelse(rhc$sex_Female == 1, "Female", "Male")
  for (v in c("edu", "das2d3pc", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1",
              "pafi1", "alb1", "hema1", "bili1", "crea1", "sod1", "pot1", "paco21", "ph1",
              "wtkilo1")) cohort[[v]] <- rhc[[v]]
  cohort$dnr1 <- yes_no("dnr1_Yes")
  cohort$ninsclas <- from_dummies(c(ninsclas_Medicaid = "Medicaid", ninsclas_Medicare = "Medicare",
                                    ninsclas_Medicare_and_Medicaid = "Medicare & Medicaid",
                                    ninsclas_No_insurance = "No insurance",
                                    ninsclas_Private_and_Medicare = "Private & Medicare"),
                                  "Private")
  for (v in c("resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps", "trauma",
              "ortho")) cohort[[v]] <- yes_no(paste0(v, "_Yes"))
  cohort$race <- ifelse(rhc$race_black == 1, "black", ifelse(rhc$race_other == 1, "other", "white"))
  cohort$income <- from_dummies(c(income1 = "$11-$25k", income2 = "$25-$50k", income3 = "> $50k"),
                                "Under $11k")
  for (v in names(cohort)) if (is.character(cohort[[v]])) cohort[[v]] <- factor(cohort[[v]])
  cohort
}

# ---- death within 30 days from the raw public file ------------------------------------------
sha256_of <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(file = path, algo = "sha256"))
  sub(" .*", "", system2("sha256sum", shQuote(path), stdout = TRUE))
}
read_death_within_30_days <- function(cohort) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  utils::download.file(rhc_url, path, mode = "wb", quiet = TRUE)
  checksum <- sha256_of(path)
  if (!identical(checksum, rhc_sha256)) {
    stop("rhc.csv has SHA-256 ", checksum, ", not the recorded ", rhc_sha256, ".")
  }
  raw <- utils::read.csv(path, stringsAsFactors = FALSE)
  # The rows are in the same order as ATbounds::RHC; check before taking the outcome.
  stopifnot(nrow(raw) == nrow(cohort),
            identical(as.numeric(raw$swang1 == "RHC"), as.numeric(cohort$swang1)),
            isTRUE(all.equal(raw$age, cohort$age)), isTRUE(all.equal(raw$aps1, cohort$aps1)))
  as.numeric(raw$dth30 == "Yes")
}

cohort <- rebuild_rhc_covariates()
covariates <- setdiff(names(cohort), "swang1")
formula <- stats::as.formula(paste("swang1 ~", paste(covariates, collapse = " + ")))

# ---- the diagnostic (cached) --------------------------------------------------------------------
if (file.exists(cache_path)) {
  sensitivity <- readRDS(cache_path)
} else {
  sensitivity <- propensity_sensitivity(formula, data = cohort,
                                        algorithms = c("logistic_regression", "bart",
                                                       "grf_probability"),
                                        replace = FALSE, caliper = 0.2, cores = 1, seed = 1)
  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(sensitivity, cache_path)
}
print(sensitivity)
used <- sensitivity$summary$models_used
treated_rows <- which(sensitivity$patients$treated == 1)
rank_disagreement <- mean(sensitivity$patients$rank_disagreement[treated_rows])
stopifnot(length(used) == 4, round(rank_disagreement, 3) == 0.047)

# ---- odds ratio of death within 30 days in each model's matched pairs --------------------------
# 1:1 matched pairs: the conditional odds ratio is the ratio of the two kinds of discordant pairs
# (treated patient died and untreated partner survived, b; the reverse, c), with the Wald interval on the log scale.
dth30 <- read_death_within_30_days(cohort)
match_set <- sensitivity$components$match_set
pair_odds_ratio <- function(id) {
  pairs <- match_set$matches[[id]]
  died_treated <- dth30[match(pairs$treated_id, match_set$row_ids)]
  died_partner <- dth30[match(pairs$control_id, match_set$row_ids)]
  b <- sum(died_treated == 1 & died_partner == 0)
  c <- sum(died_treated == 0 & died_partner == 1)
  se <- sqrt(1 / b + 1 / c)
  data.frame(model_id = id, matched_pairs = nrow(pairs), odds_ratio = b / c,
             conf_low = exp(log(b / c) - stats::qnorm(0.975) * se),
             conf_high = exp(log(b / c) + stats::qnorm(0.975) * se))
}
models <- sensitivity$models[match(match_set$model_ids, sensitivity$models$model_id),
                             c("model_id", "algorithm", "logloss", "c_statistic", "smd_max_main",
                               "smd_max_extended", "treated_kept", "passed_logloss_gate",
                               "passed_quality_gate")]
models <- cbind(models, do.call(rbind, lapply(models$model_id, pair_odds_ratio))[, -1])
rownames(models) <- NULL
utils::write.csv(models, file.path(figure_dir, "rhc-models.csv"), row.names = FALSE)
print(models)

# ---- figures ------------------------------------------------------------------------------
okabe <- c(blue = "#0072B2", vermillion = "#D55E00", green = "#009E73", orange = "#E69F00",
           sky = "#56B4E9", grey = "#999999", black = "#000000")
model_label <- function(id) {
  algorithm <- sub("_[0-9]+$", "", id)
  ifelse(algorithm == "logistic_regression", "Logistic regression",
         paste0(ifelse(algorithm == "bart", "BART", "GRF"), ", seed ", sub("^.*_", "", id)))
}
theme_example <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title.position = "plot",
          plot.title = element_text(face = "bold"), legend.position = "bottom")
}
save_figure <- function(plot, name, height) {
  ggsave(file.path(figure_dir, name), plot, width = 8, height = height, dpi = 100, bg = "white")
}
model_set <- sensitivity$components$model_set
row_of <- function(ids) match(ids, match_set$row_ids)

# Figure 1: one treated patient and its matched partner under each kept model. The patient is
# one matched by all four models with the most distinct partners and, among those, the rank
# disagreement closest to the study mean (a typical patient, not an extreme one).
partner <- lapply(stats::setNames(used, used), function(id) {
  pairs <- match_set$matches[[id]]
  stats::setNames(as.character(pairs$control_id), as.character(pairs$treated_id))
})
matched_by_all <- Reduce(intersect, lapply(partner, names))
distinct <- vapply(matched_by_all, function(t) length(unique(vapply(partner, `[[`, "", t))), 1L)
patient_d <- sensitivity$patients$rank_disagreement[match(matched_by_all, sensitivity$patients$id)]
candidates <- which(distinct == max(distinct))
chosen <- matched_by_all[candidates[which.min(abs(patient_d[candidates] - rank_disagreement))]]
partners <- vapply(partner, `[[`, "", chosen)
partner_letter <- stats::setNames(paste("Untreated patient", LETTERS[seq_along(unique(partners))]),
                                  unique(partners))
figure_1_data <- data.frame(
  model = factor(model_label(used), levels = rev(model_label(used))),
  treated_rank = 100 * model_set$rank[row_of(chosen), used],
  partner_rank = 100 * model_set$rank[cbind(row_of(partners),
                                            match(used, colnames(model_set$rank)))],
  partner = unname(partner_letter[partners]))
figure_1_data$y <- rev(seq_len(nrow(figure_1_data)))
figure_1_data$model_text <- sprintf("%s\ntreated patient at rank %.1f", as.character(figure_1_data$model),
                                    figure_1_data$treated_rank)
figure_1_data$partner_text <- sprintf("%s\nrank %.1f", figure_1_data$partner,
                                      figure_1_data$partner_rank)
centre <- mean(figure_1_data$y)
partner_colours <- stats::setNames(unname(okabe[c("vermillion", "blue", "green", "orange")])[
  seq_along(partner_letter)], unname(partner_letter))
figure_1 <- ggplot(figure_1_data) +
  geom_segment(aes(x = 0.18, xend = 0.82, y = centre, yend = y), colour = okabe[["grey"]],
               linewidth = 0.8) +
  geom_segment(aes(x = 1.18, xend = 1.78, y = y, yend = y, colour = partner), linewidth = 1,
               arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")) +
  annotate("label", x = 0, y = centre, label = "Treated patient\n(catheterized)", size = 4.2,
           fontface = "bold", label.padding = grid::unit(0.5, "lines"), fill = "white") +
  geom_label(aes(x = 1, y = y, label = model_text), size = 3.5, fill = "grey96",
             linewidth = 0) +
  geom_label(aes(x = 2.05, y = y, label = partner_text, colour = partner), size = 3.6,
             fontface = "bold", fill = "white", show.legend = FALSE) +
  scale_colour_manual(values = partner_colours, guide = "none") +
  scale_x_continuous(limits = c(-0.25, 2.35)) +
  scale_y_continuous(limits = c(0.5, nrow(figure_1_data) + 0.5)) +
  labs(title = paste0("One treated patient, four kept models, ",
                      c("one", "two", "three", "four")[max(distinct)],
                      " different untreated partners"),
       subtitle = paste0("Each model ranks the patient a little differently (percentile rank of ",
                         "the propensity score, 0 to 100)\nand chooses the nearest untreated ",
                         "patient on its own scale.")) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.margin = margin(8, 8, 8, 8))
save_figure(figure_1, "rhc-partners.png", 4)

# Study-level numbers quoted in the vignette.
table_value <- function(name) sensitivity$summary$table$with_gate[
  sensitivity$summary$table$quantity == name]
study <- data.frame(
  quantity = c("patients", "treated", "models", "models_kept", "rank_disagreement",
               "rank_disagreement_between", "rank_disagreement_within", "pair_instability",
               "bootstrap_rank_instability", "smd_max_main_range", "treated_kept_min",
               "treated_kept_max", "treated_matched_by_all_kept",
               "share_with_more_than_one_partner", "example_distinct_partners"),
  value = c(nrow(cohort), length(treated_rows), length(match_set$model_ids), length(used),
            rank_disagreement,
            mean(sensitivity$patients$rank_disagreement_between[treated_rows]),
            mean(sensitivity$patients$rank_disagreement_within[treated_rows]),
            table_value("pair_instability_models_mean_treated"),
            table_value("bootstrap_rank_instability_treated"),
            table_value("smd_max_main_range"),
            min(models$treated_kept), max(models$treated_kept), length(matched_by_all),
            mean(distinct > 1), max(distinct)))
utils::write.csv(study, file.path(figure_dir, "rhc-study.csv"), row.names = FALSE)
print(study)

# Figure 2: the quality gate.
gate <- models
gate$label <- factor(model_label(gate$model_id), levels = rev(model_label(gate$model_id)))
gate$status <- ifelse(gate$passed_quality_gate, "Kept", "Excluded")
limits <- c(min(gate$logloss) * (1 + sensitivity$settings$quality_tolerance),
            sensitivity$settings$balance_limit)
panels <- c("Fit: cross-fitted log-loss", "Balance: largest |SMD| after matching")
figure_2_data <- rbind(
  data.frame(label = gate$label, status = gate$status, value = gate$logloss, panel = panels[1],
             limit = limits[1]),
  data.frame(label = gate$label, status = gate$status, value = gate$smd_max_main,
             panel = panels[2], limit = limits[2]))
figure_2_data$panel <- factor(figure_2_data$panel, levels = panels)
figure_2 <- ggplot(figure_2_data, aes(value, label, colour = status, shape = status)) +
  geom_vline(aes(xintercept = limit), linetype = "dashed", colour = okabe[["grey"]]) +
  geom_point(size = 3.5, stroke = 1.3) +
  facet_wrap(~panel, scales = "free_x") +
  scale_colour_manual(values = c(Kept = okabe[["blue"]], Excluded = okabe[["vermillion"]]),
                      name = NULL) +
  scale_shape_manual(values = c(Kept = 16, Excluded = 4), name = NULL) +
  labs(x = "Lower is better. Dashed line: the limit of the gate",
       y = NULL, title = "The quality gate: four of seven models kept",
       subtitle = "A model is kept when it is at or left of the dashed line in both panels") +
  theme_example()
save_figure(figure_2, "rhc-gate.png", 3.8)

# Figure 3: how far the treated patients' ranks move between the kept models.
bart_used <- setdiff(used, "logistic_regression_1")
figure_3_data <- do.call(rbind, lapply(bart_used, function(id) {
  data.frame(logistic = 100 * model_set$rank[treated_rows, "logistic_regression_1"],
             other = 100 * model_set$rank[treated_rows, id], model = model_label(id))
}))
figure_3 <- ggplot(figure_3_data, aes(logistic, other)) +
  geom_abline(slope = 1, intercept = 0, colour = okabe[["grey"]], linetype = "dashed") +
  geom_point(alpha = 0.25, size = 0.7, colour = okabe[["blue"]]) +
  facet_wrap(~model) +
  coord_equal() +
  labs(x = "Percentile rank under logistic regression", y = "Percentile rank under BART",
       title = "The same catheterized patients, ranked by different models",
       subtitle = paste0("Each dot is one of the ", format(length(treated_rows), big.mark = ","),
                         " catheterized patients; on the dashed line both models agree.\n",
                         "Mean rank disagreement over the four kept models: ",
                         sprintf("%.3f", rank_disagreement), ".")) +
  theme_example()
save_figure(figure_3, "rhc-ranks.png", 4.2)

# Figure 4: odds ratio of death within 30 days in each kept model's matched pairs.
kept <- models[models$model_id %in% used, ]
kept$label <- factor(model_label(kept$model_id), levels = rev(model_label(kept$model_id)))
kept$balance <- ifelse(round(kept$smd_max_main, 2) <= 0.10,
                       "Largest |SMD| after matching at most 0.10 (two decimals)",
                       "Largest |SMD| after matching slightly above 0.10")
figure_4 <- ggplot(kept, aes(odds_ratio, label, colour = balance, shape = balance)) +
  geom_vline(xintercept = 1, colour = okabe[["grey"]]) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2, orientation = "y") +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%.2f", odds_ratio)), vjust = -1.1, size = 3.6,
            show.legend = FALSE) +
  scale_colour_manual(values = unname(okabe[c("blue", "orange")]), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_log10(breaks = c(1, 1.2, 1.4, 1.6), limits = c(0.95, 1.65)) +
  labs(x = "Odds ratio of death within 30 days, catheterized vs matched untreated (95 % CI)",
       y = NULL, title = "Effect estimates across the four kept models",
       subtitle = "Grey line: odds ratio 1, no difference") +
  theme_example()
save_figure(figure_4, "rhc-effects.png", 3.6)

print(file.info(file.path(figure_dir, c("rhc-partners.png", "rhc-gate.png",
                                        "rhc-ranks.png", "rhc-effects.png")))["size"])
