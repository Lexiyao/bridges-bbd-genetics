# =============================================================================
# Dissertation Analysis — Objective 2: Hormonal and Reproductive Factors
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil
#
# DATA GOVERNANCE
#   Data accessed under BCAC/BRIDGES approved secondary use agreement.
#   All analysis conducted on non-identifiable summary-level outputs only.
#
# BBD DEFINITION (corrected):
#   BBD cases    : status == 0 & ethnicityClass == 1 & BBD_history == 1 & !is.na(ageInt)
#   BBD controls : status == 0 & ethnicityClass == 1 & BBD_history == 0 & !is.na(ageInt)
#   (Previous first-attempt script incorrectly used status == 3 for BBD — WRONG for
#    this dataset. status == 3 in BRIDGES v17 has 89% missing BBD_history and does
#    not correspond to the BBD phenotype.)
#
# DCIS controls / LCIS controls: all status==0 & ethnicityClass==1 & !is.na(ageInt)
#   consistent with missense_analysis_corrected.R
#
# MODEL: standard logistic regression (glm, binomial)
#   Hormonal/reproductive exposures are common (not rare), so standard glm is
#   appropriate. Firth's regression is reserved for rare-variant analyses (Obj 1).
#
# COVARIATES: ageInt (continuous) + study (categorical)
# FDR: Benjamini-Hochberg across all 15 tests (5 exposures × 3 outcomes)
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c("tidyverse", "broom", "writexl", "ggplot2", "scales")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(tidyverse)
library(broom)
library(writexl)
library(ggplot2)
library(scales)

# ── 1. FILE PATH ──────────────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"

# Ensure output directories exist
dir.create("outputs/tables",       recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures",      recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/session_info", recursive = TRUE, showWarnings = FALSE)

# Working-directory guard: scripts use paths relative to the project root.
if (!file.exists(PHENO_FILE))
  stop("Cannot find '", PHENO_FILE, "'. Run this script from the project root ",
       "(the directory containing the BRIDGES data files).", call. = FALSE)

# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
message("Loading phenotype file...")
pheno <- read_delim(
  PHENO_FILE, delim = "\t", show_col_types = FALSE,
  na = c("", "NA", "888", "777", "999")
)
message(sprintf("Loaded: %d rows", nrow(pheno)))

# ── 2b. BINARY EXPOSURE VALUE CHECK ───────────────────────────────────────────
# Print observed values for menoStat and HRTEver at runtime.
# menoStat is coded 1/2 (NOT 0/1); HRTEver is coded 0/1.
# Expected values: menoStat in {1, 2, NA}; HRTEver in {0, 1, NA}.
# If unexpected values appear, verify against BCAC Extended Data Dictionary v4.0
# and add any unlisted codes to the na= list in read_delim() above.
message("\nValue check — menoStat (expected: 1=pre/peri, 2=post, NA):")
print(table(pheno$menoStat, useNA = "always"))
message("Value check — HRTEver (expected: 0=never, 1=ever, NA):")
print(table(pheno$HRTEver,  useNA = "always"))

# ── 3. EXPOSURES ──────────────────────────────────────────────────────────────
# parity      : number of full-term pregnancies (continuous)
# ageFFTP     : age at first full-term pregnancy (continuous)
#               NOTE: ageFFTP is structurally missing for nulliparous women.
#               drop_na() inside run_glm() therefore silently restricts the
#               ageFFTP models to parous women only. The effective N for these
#               models will be lower than for other exposures and is reported
#               via n_complete in the output. This restriction is noted in
#               dissertation Results (Section 4.5).
# ageMenarche : age at menarche (continuous)
# menoStat    : menopausal status (1=pre/peri-menopausal, 2=post-menopausal)
#               NOTE: coded 1/2 in BRIDGES v17 (NOT 0/1). Confirmed against
#               BCAC Extended Data Dictionary v4.0. Modelled as continuous;
#               OR represents the contrast post (2) vs pre (1), equivalent
#               in magnitude to a 0/1 binary coding. See dictionary: MenoStat.
# HRTEver     : ever HRT use (0=never, 1=ever)
EXPOSURES <- c("parity", "ageFFTP", "ageMenarche", "menoStat", "HRTEver")

# ── 4. BUILD COHORTS ──────────────────────────────────────────────────────────
# European, cancer-free, non-missing age
eur_status0 <- pheno %>%
  filter(ethnicityClass == 1, status == 0, !is.na(ageInt))

message(sprintf("European status-0 pool: %d", nrow(eur_status0)))

# BBD cohort (CORRECTED definition)
bbd_cases <- eur_status0 %>%
  filter(BBD_history == 1) %>%
  mutate(outcome = 1L)

bbd_ctrls <- eur_status0 %>%
  filter(BBD_history == 0) %>%
  mutate(outcome = 0L)

dat_bbd <- bind_rows(bbd_cases, bbd_ctrls)

message(sprintf("BBD cases   : %d  [expected 1,765]", nrow(bbd_cases)))
message(sprintf("BBD controls: %d  [expected 12,715]", nrow(bbd_ctrls)))

# DCIS cohort
dcis_cases <- pheno %>%
  filter(ethnicityClass == 1, status == 2,
         MorphologygroupIndex_corr == "Ductal", !is.na(ageInt)) %>%
  mutate(outcome = 1L)

dcis_ctrls <- eur_status0 %>% mutate(outcome = 0L)
dat_dcis <- bind_rows(dcis_ctrls, dcis_cases)

message(sprintf("DCIS cases   : %d  [expected 1,663]", nrow(dcis_cases)))
message(sprintf("DCIS controls: %d  [expected 34,424]", nrow(dcis_ctrls)))

# LCIS cohort
lcis_cases <- pheno %>%
  filter(ethnicityClass == 1, status == 2,
         MorphologygroupIndex_corr == "Lobular", !is.na(ageInt)) %>%
  mutate(outcome = 1L)

lcis_ctrls <- eur_status0 %>% mutate(outcome = 0L)
dat_lcis <- bind_rows(lcis_ctrls, lcis_cases)

message(sprintf("LCIS cases   : %d  [expected 139]", nrow(lcis_cases)))
message(sprintf("LCIS controls: %d  [expected 34,424]", nrow(lcis_ctrls)))

# ── 5. REGRESSION FUNCTION ────────────────────────────────────────────────────
run_glm <- function(data, exposure, outcome_label) {
  dat_sub <- data %>%
    select(outcome, all_of(exposure), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()

  n_cases    <- sum(dat_sub$outcome == 1)
  n_complete <- nrow(dat_sub)

  formula_str <- paste("outcome ~", exposure, "+ ageInt + study")

  fit <- tryCatch(
    glm(as.formula(formula_str), data = dat_sub, family = binomial()),
    error = function(e) {
      message(sprintf("  glm failed for %s / %s: %s", exposure, outcome_label, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) {
    return(tibble(exposure, outcome = outcome_label,
                  n_cases, n_complete,
                  estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                  p.value = NA_real_))
  }

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure) %>%
    mutate(exposure = exposure, outcome = outcome_label,
           n_cases = n_cases, n_complete = n_complete) %>%
    select(exposure, outcome, n_cases, n_complete,
           estimate, conf.low, conf.high, p.value)
}

# ── 6. RUN ALL 15 MODELS ──────────────────────────────────────────────────────
message("\nRunning 15 logistic regression models...")

datasets <- list(BBD = dat_bbd, DCIS = dat_dcis, LCIS = dat_lcis)

results_obj2 <- map_dfr(names(datasets), function(out_label) {
  map_dfr(EXPOSURES, function(exp) {
    message(sprintf("  %s / %s", exp, out_label))
    run_glm(datasets[[out_label]], exp, out_label)
  })
})

# BH-FDR across all 15 tests
results_obj2 <- results_obj2 %>%
  mutate(p.fdr = p.adjust(p.value, method = "BH"),
         fdr_sig = p.fdr < 0.05)

message("\n── Results ──")
print(results_obj2, n = Inf)

# ── 7. FORMAT TABLE ───────────────────────────────────────────────────────────
exposure_labels <- c(
  parity      = "Parity (per additional birth)",
  ageFFTP     = "Age at first full-term pregnancy (per year)",
  ageMenarche = "Age at menarche (per year)",
  menoStat    = "Post-menopausal vs pre-menopausal",
  HRTEver     = "Ever HRT use vs never"
)

fmt_obj2 <- results_obj2 %>%
  mutate(
    Exposure   = exposure_labels[exposure],
    Outcome    = outcome,
    `N cases`  = n_cases,
    `N (complete)` = n_complete,
    `OR (95% CI)` = case_when(
      is.na(estimate) ~ "—",
      TRUE ~ sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)
    ),
    `p-value` = case_when(
      is.na(p.value)  ~ "—",
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    ),
    `p-FDR (BH)` = case_when(
      is.na(p.fdr)  ~ "—",
      p.fdr < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.fdr)
    ),
    `FDR sig.` = if_else(fdr_sig, "Yes*", "No", missing = "—")
  ) %>%
  select(Exposure, Outcome, `N cases`, `N (complete)`,
         `OR (95% CI)`, `p-value`, `p-FDR (BH)`, `FDR sig.`)

# ── 8. EXPORT EXCEL ───────────────────────────────────────────────────────────
out_path <- "outputs/tables/objective2_hormonal_corrected.xlsx"
write_xlsx(
  list(
    "Hormonal_results" = fmt_obj2,
    "Raw"              = results_obj2
  ),
  path = out_path
)
message(sprintf("\nExported: %s", out_path))

# ── 9. FOREST PLOT ────────────────────────────────────────────────────────────
plot_df <- results_obj2 %>%
  filter(!is.na(estimate), conf.high < 15) %>%
  mutate(
    outcome_f  = factor(outcome, levels = c("BBD", "DCIS", "LCIS")),
    exposure_f = factor(exposure,
                        levels = rev(c("parity","ageFFTP","ageMenarche",
                                       "menoStat","HRTEver")),
                        labels = rev(c("Parity","Age at FFTP",
                                       "Age at menarche",
                                       "Post-menopausal","Ever HRT"))),
    sig_group  = if_else(fdr_sig, "FDR < 0.05", "Not significant")
  )

source("scripts/_thesis_theme.R")   # THESIS_PAL, theme_thesis(), save_fig()
p_obj2 <- ggplot(plot_df, aes(x = estimate, y = exposure_f, colour = sig_group)) +
  geom_vline(xintercept = 1, linetype = "longdash",
             colour = unname(THESIS_PAL["rule"]), linewidth = 0.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.5, shape = 18) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2),
                labels = c("0.25", "0.5", "1", "2")) +
  scale_colour_manual(values = c("FDR < 0.05"     = unname(THESIS_PAL["signal"]),
                                 "Not significant" = unname(THESIS_PAL["baseline"]))) +
  facet_wrap(~ outcome_f, nrow = 1) +
  labs(
    title   = "Hormonal and reproductive factor associations with BBD, DCIS and LCIS",
    x       = "Odds Ratio (95% CI, log scale)",
    y       = NULL,
    colour  = NULL,
    caption = "Adjusted for age at interview and study. BH-FDR across 15 tests."
  ) +
  theme_thesis(base_size = 13, grid = "y") +
  theme(strip.background = element_rect(fill = "#EAF2FF", colour = NA),
        strip.text       = element_text(face = "bold"))

save_fig(p_obj2, "outputs/figures/forest_hormonal_corrected.png", width = 12, height = 5)

# ── 10. SESSION INFO ──────────────────────────────────────────────────────────
sink("outputs/session_info/session_info_obj2_hormonal.txt")
cat("Objective 2 Hormonal Analysis — Session Info\n")
cat(format(Sys.time()), "\n\n")
sessionInfo()
sink()
message("Session info saved.")
