# =============================================================================
# Missingness Audit — Objective 2 hormonal/reproductive complete-case analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE (addresses supervisor point #2: "the hormonal objective uses
#          complete-case analysis — check this")
#   objective2_hormonal_analysis.R fits each exposure with drop_na(), i.e.
#   complete-case. Complete-case analysis is UNBIASED if data are missing
#   completely at random (MCAR) or missingness is unrelated to the outcome
#   conditional on covariates. It is BIASED if missingness differs by
#   case/control status. This script quantifies, for every exposure × outcome:
#     (i)  the % missing overall and separately in cases vs controls
#     (ii) a test for DIFFERENTIAL missingness by case status (chi-square)
#   so the complete-case assumption can be defended or flagged per exposure.
#
#   ageFFTP is STRUCTURALLY missing for nulliparous women (not informative
#   missingness) and is flagged separately.
#
# DATA GOVERNANCE: BCAC/BRIDGES approved secondary use; non-disclosive output.
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(writexl)
})

PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

message("Loading phenotype file...")
pheno <- suppressWarnings(read_delim(
  PHENO_FILE, delim = "\t", show_col_types = FALSE,
  na = c("", "NA", "888", "777", "999")))

EXPOSURES <- c("parity", "ageFFTP", "ageMenarche", "menoStat", "HRTEver")

eur_status0 <- pheno %>% filter(ethnicityClass == 1, status == 0, !is.na(ageInt))

# Cohorts identical to objective2_hormonal_analysis.R
bbd  <- bind_rows(
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

# For one cohort × exposure: missingness in cases vs controls + differential test
audit_one <- function(data, exposure, outcome_label) {
  d <- data %>% mutate(miss = is.na(.data[[exposure]]))
  n_case <- sum(d$outcome == 1); n_ctrl <- sum(d$outcome == 0)
  miss_case <- sum(d$miss & d$outcome == 1)
  miss_ctrl <- sum(d$miss & d$outcome == 0)

  # chi-square test: is missingness associated with case status?
  tab <- table(d$outcome, d$miss)
  p_diff <- tryCatch(suppressWarnings(chisq.test(tab)$p.value),
                     error = function(e) NA_real_)

  tibble(
    Outcome = outcome_label, Exposure = exposure,
    `% miss (cases)` = 100 * miss_case / n_case,
    `% miss (ctrls)` = 100 * miss_ctrl / n_ctrl,
    `% miss (all)`   = 100 * (miss_case + miss_ctrl) / (n_case + n_ctrl),
    `N analysed (cases)` = n_case - miss_case,
    `N dropped (cases)`  = miss_case,
    `p (differential)`   = p_diff
  )
}

audit <- map_dfr(names(cohorts), function(oc) {
  map_dfr(EXPOSURES, function(ex) audit_one(cohorts[[oc]], ex, oc))
})

audit <- audit %>%
  mutate(
    structural = if_else(Exposure == "ageFFTP",
                         "YES — missing for nulliparous (not informative)", ""),
    differential_flag = case_when(
      Exposure == "ageFFTP"        ~ "n/a (structural)",
      `p (differential)` < 0.001   ~ "STRONG differential missingness",
      `p (differential)` < 0.05    ~ "differential missingness",
      TRUE                          ~ "no evidence of differential missingness"
    )
  )

options(pillar.sigfig = 3)
cat("\n", strrep("=", 92), "\n",
    "MISSINGNESS AUDIT — hormonal/reproductive exposures (complete-case diagnostic)\n",
    strrep("=", 92), "\n", sep = "")
audit %>%
  mutate(across(starts_with("% miss"), ~ sprintf("%.1f", .x)),
         `p (differential)` = if_else(`p (differential)` < 0.001, "<0.001",
                                      sprintf("%.3f", `p (differential)`))) %>%
  select(Outcome, Exposure, `% miss (cases)`, `% miss (ctrls)`,
         `p (differential)`, differential_flag) %>%
  as.data.frame() %>% print(right = FALSE)

write_xlsx(list("Missingness_audit" = as.data.frame(audit)),
           "outputs/tables/hormonal_missingness_audit.xlsx")
message("\nExported: outputs/tables/hormonal_missingness_audit.xlsx")
