# =============================================================================
# Age-Model Sensitivity — menopausal status & HRT (Objective 2)
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE (addresses supervisor point #7: "the age model affects not just
#          menopause but also HRT")
#   The primary Objective-2 models adjust for age LINEARLY (ageInt). Both
#   menopausal status and ever-HRT use are strongly age-dependent, so a
#   mis-specified (too rigid) age adjustment could leave residual confounding
#   that masquerades as a menoStat or HRT effect. This script refits the
#   menoStat and HRTEver models for all three outcomes under three age
#   specifications and checks whether the exposure OR is stable:
#     (L) linear age            ageInt              [= PRIMARY, as reported]
#     (Q) quadratic age         ageInt + ageInt^2
#     (S) natural spline age    ns(ageInt, df = 3)
#
#   A stable OR across (L)-(Q)-(S) => the association is not an artifact of
#   linear-age mis-specification. A large shift => residual age confounding.
#
# Cohorts & complete-case handling identical to objective2_hormonal_analysis.R.
# DATA GOVERNANCE: BCAC/BRIDGES approved secondary use; non-disclosive output.
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(broom)
  library(splines)
  library(writexl)
})

PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

pheno <- suppressWarnings(read_delim(
  PHENO_FILE, delim = "\t", show_col_types = FALSE,
  na = c("", "NA", "888", "777", "999")))

eur_status0 <- pheno %>% filter(ethnicityClass == 1, status == 0, !is.na(ageInt))

bbd <- bind_rows(
  eur_status0 %>% filter(BBD_history == 1) %>% mutate(outcome = 1L),
  eur_status0 %>% filter(BBD_history == 0) %>% mutate(outcome = 0L))
dcis <- bind_rows(
  eur_status0 %>% mutate(outcome = 0L),
  pheno %>% filter(ethnicityClass == 1, status == 2,
                   MorphologygroupIndex_corr == "Ductal", !is.na(ageInt)) %>%
    mutate(outcome = 1L))
lcis <- bind_rows(
  eur_status0 %>% mutate(outcome = 0L),
  pheno %>% filter(ethnicityClass == 1, status == 2,
                   MorphologygroupIndex_corr == "Lobular", !is.na(ageInt)) %>%
    mutate(outcome = 1L))
cohorts <- list(BBD = bbd, DCIS = dcis, LCIS = lcis)

AGE_TERMS <- c(L = "ageInt",
               Q = "ageInt + I(ageInt^2)",
               S = "ns(ageInt, df = 3)")

fit_age <- function(data, exposure, outcome_label, age_key) {
  d <- data %>%
    select(outcome, all_of(exposure), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()
  form <- as.formula(paste("outcome ~", exposure, "+", AGE_TERMS[age_key], "+ study"))
  fit  <- glm(form, data = d, family = binomial())
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure) %>%
    transmute(outcome = outcome_label, exposure, age_model = age_key,
              n_complete = nrow(d),
              OR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value)
}

EXPOSURES <- c("menoStat", "HRTEver")
grid <- expand_grid(outcome = names(cohorts), exposure = EXPOSURES,
                    age_key = names(AGE_TERMS))

res <- pmap_dfr(grid, function(outcome, exposure, age_key)
  fit_age(cohorts[[outcome]], exposure, outcome, age_key))

wide <- res %>%
  mutate(cell = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)) %>%
  select(outcome, exposure, age_model, cell, p) %>%
  pivot_wider(names_from = age_model, values_from = c(cell, p)) %>%
  transmute(
    Outcome = outcome, Exposure = exposure,
    `OR linear (PRIMARY)` = cell_L, `p linear` = sprintf("%.3f", p_L),
    `OR quadratic`        = cell_Q,
    `OR spline df3`       = cell_S,
    # max absolute fractional change in OR vs the linear estimate
    note = "")

# attach a robustness flag computed from the numeric ORs
flag <- res %>%
  group_by(outcome, exposure) %>%
  summarise(OR_L = OR[age_model == "L"],
            max_rel_shift = max(abs(OR[age_model != "L"] - OR[age_model == "L"]) /
                                OR[age_model == "L"]),
            .groups = "drop") %>%
  mutate(stability = if_else(max_rel_shift < 0.10,
                             sprintf("STABLE (<10%% shift: %.1f%%)", 100 * max_rel_shift),
                             sprintf("SHIFT %.1f%% — check residual age confounding",
                                     100 * max_rel_shift)))

wide <- wide %>% left_join(flag %>% select(outcome, exposure, stability),
                           by = c("Outcome" = "outcome", "Exposure" = "exposure")) %>%
  select(-note)

cat("\n", strrep("=", 96), "\n",
    "AGE-MODEL SENSITIVITY — menoStat & HRTEver under linear / quadratic / spline age\n",
    strrep("=", 96), "\n", sep = "")
as.data.frame(wide) %>% print(right = FALSE)

write_xlsx(list("Age_sensitivity" = as.data.frame(wide),
                "Raw" = as.data.frame(res)),
           "outputs/tables/age_model_sensitivity.xlsx")
message("\nExported: outputs/tables/age_model_sensitivity.xlsx")
