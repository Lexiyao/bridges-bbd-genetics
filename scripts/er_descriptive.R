# =============================================================================
# Descriptive ER-status distribution of carriers among DCIS cases
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — University of Cambridge MPhil
#
# PURPOSE
#   ER status (ER_statusIndex) is available for ~41% of DCIS cases (678/1663).
#   A formal
#   ER-stratified case-control regression is not pursued because the ER-negative
#   subgroup and per-gene carrier counts are too sparse (would reproduce the
#   penalised-estimate artifacts seen for LCIS). Instead, this descriptive
#   analysis asks, among DCIS cases with KNOWN ER status, whether carriers of
#   the FDR-significant genes are more likely to be ER-positive than non-carriers
#   — a low-risk check that bears on the ER-status interpretation in the
#   Discussion. Fisher's exact test is reported but is exploratory and
#   underpowered given small carrier numbers and the high baseline ER-positivity
#   of DCIS.
#
# OUTPUT
#   outputs/tables/er_descriptive_carriers.xlsx
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(writexl)})
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# NB: na must include 888/777/999 so that cases with unknown age (ageInt==888)
# are dropped by the !is.na(ageInt) filter, matching the main analysis cohort
# (DCIS n = 1,663). Reading with na=c("","NA") alone retains 118 unknown-age
# cases and inflates the DCIS denominator to 1,781, which is inconsistent with
# every regression in the study. ER_statusIndex==888 (missing) is likewise
# nulled here and correctly excluded by the %in% c(0,1) filter below.
pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                    show_col_types = FALSE, na = c("", "NA", "888", "777", "999"))
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv", show_col_types = FALSE, na = c("","NA"))
miss  <- read_csv("concept_807_zhang_bridges_missense.csv",  show_col_types = FALSE, na = c("","NA"))

# DCIS cases with KNOWN ER status (ER_statusIndex 1=positive, 0=negative; 888=missing)
dcis <- pheno %>%
  filter(status == 2, ethnicityClass == 1, !is.na(ageInt),
         MorphologygroupIndex_corr == "Ductal",
         ER_statusIndex %in% c(0, 1)) %>%
  mutate(ER_pos = as.integer(ER_statusIndex == 1)) %>%
  left_join(trunc %>% select(BRIDGES_ID, ATM_truncating, BRCA2_truncating, CHEK2_truncating), by = "BRIDGES_ID") %>%
  left_join(miss  %>% select(BRIDGES_ID, CHEK2_CADD.phred.01, TP53_CADD.phred.01), by = "BRIDGES_ID")

baseline_erpos <- mean(dcis$ER_pos)   # overall ER+ proportion among DCIS with known ER
cat(sprintf("DCIS with known ER status: n = %d (%.0f%% ER-positive)\n\n",
            nrow(dcis), 100 * baseline_erpos))

genes <- tribble(
  ~label,              ~col,
  "ATM (truncating)",   "ATM_truncating",
  "BRCA2 (truncating)", "BRCA2_truncating",
  "CHEK2 (truncating)", "CHEK2_truncating",
  "CHEK2 (missense)",   "CHEK2_CADD.phred.01",
  "TP53 (missense)",    "TP53_CADD.phred.01"
)

res <- pmap_dfr(genes, function(label, col) {
  carr <- as.integer(replace_na(as.numeric(dcis[[col]]), 0) >= 1)
  # 2x2: carrier (yes/no) x ER (pos/neg)
  cp <- sum(carr == 1 & dcis$ER_pos == 1); cn <- sum(carr == 1 & dcis$ER_pos == 0)
  np <- sum(carr == 0 & dcis$ER_pos == 1); nn <- sum(carr == 0 & dcis$ER_pos == 0)
  n_carr <- cp + cn
  fish <- if (n_carr > 0) fisher.test(matrix(c(cp, cn, np, nn), nrow = 2))$p.value else NA_real_
  tibble(Gene = label,
         `Carriers (known ER)` = n_carr,
         `Carrier ER+ n` = cp, `Carrier ER- n` = cn,
         `Carrier ER+ %` = if (n_carr > 0) round(100 * cp / n_carr) else NA,
         `Non-carrier ER+ %` = round(100 * np / (np + nn)),
         `Fisher p` = round(fish, 3))
})

cat("Carrier ER-status distribution among DCIS cases (known ER only):\n")
print(as.data.frame(res))
cat(sprintf("\nBaseline: %.0f%% of DCIS (known ER) are ER-positive.\n", 100 * baseline_erpos))
cat("Interpretation: carriers are predominantly ER-positive, consistent with\n")
cat("the ER-status hypothesis; numbers are small and Fisher tests exploratory.\n")

write_xlsx(list(er_carrier_distribution = res,
                meta = tibble(metric = c("DCIS with known ER", "Baseline ER+ %"),
                              value = c(nrow(dcis), round(100 * baseline_erpos)))),
           "outputs/tables/er_descriptive_carriers.xlsx")
message("\nSaved: outputs/tables/er_descriptive_carriers.xlsx")
