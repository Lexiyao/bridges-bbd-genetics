# =============================================================================
# Collider / selection-bias check for the BBD genetic analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE
#   The BBD genetic analysis defines cases as cancer-free women with a BBD
#   history (status==0 & BBD_history==1). Because the panel genes CAUSE breast
#   cancer, BBD-history women who carry a variant are more likely to have
#   progressed to cancer and thus to have left the status==0 stratum — so the
#   cancer-free BBD case group is selectively DEPLETED of carriers. Conditioning
#   the cases on cancer-free status (a descendant of the gene) is a collider,
#   biasing the BBD odds ratios toward (or below) the null.
#
#   This script (1) quantifies the depletion — carrier frequency among
#   BBD-history women who are cancer-free vs those who progressed to cancer —
#   and (2) contrasts three case definitions to show the bias is not removable
#   by simple adjustment:
#     (A) cancer-free BBD            [primary; collider-biased toward null]
#     (B) any-status BBD             [avoids collider but conflates with cancer]
#     (C) adjust for cancer status   [not estimable: cases are all status==0]
#
# OUTPUT
#   outputs/tables/bbd_collider_check.xlsx
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(writexl)})
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                    show_col_types = FALSE, na = c("", "NA", "888", "777", "999"))
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv", show_col_types = FALSE, na = c("","NA"))
GENES <- c("CHEK2", "ATM", "BRCA2", "BRCA1", "PALB2")

dat <- pheno %>% filter(ethnicityClass == 1) %>%
  left_join(trunc %>% select(BRIDGES_ID, all_of(paste0(GENES, "_truncating"))), by = "BRIDGES_ID")
carr <- function(g) as.integer(replace_na(as.numeric(dat[[paste0(g, "_truncating")]]), 0) >= 1)

# (1) Depletion: carrier frequency among BBD-history women, cancer-free vs cancer
depletion <- map_dfr(GENES, function(g) {
  c <- carr(g)
  cf <- dat$BBD_history == 1 & dat$status == 0          # cancer-free BBD (cases)
  ca <- dat$BBD_history == 1 & dat$status %in% c(1, 2)  # BBD that became cancer (excluded)
  tibble(Gene = g,
         `Cancer-free BBD: carriers` = sum(c[cf], na.rm = TRUE),
         `Cancer-free BBD: n` = sum(cf, na.rm = TRUE),
         `Cancer-free freq %` = round(100 * sum(c[cf], na.rm=TRUE) / sum(cf, na.rm=TRUE), 2),
         `Cancer BBD: carriers` = sum(c[ca], na.rm = TRUE),
         `Cancer BBD: n` = sum(ca, na.rm = TRUE),
         `Cancer freq %` = round(100 * sum(c[ca], na.rm=TRUE) / sum(ca, na.rm=TRUE), 2),
         `Fold (cancer / cancer-free)` = round((sum(c[ca],na.rm=TRUE)/sum(ca,na.rm=TRUE)) /
                                               (sum(c[cf],na.rm=TRUE)/sum(cf,na.rm=TRUE)), 1))
})

# (2) Case-definition contrast (age + study adjusted glm).
# case_mask: a logical vector selecting the cases; controls are always the
# cancer-free, BBD-free women.
or_def <- function(case_mask, g) {
  base <- dat %>% mutate(c = as.integer(replace_na(as.numeric(.data[[paste0(g, "_truncating")]]), 0) >= 1))
  cases <- base %>% filter(!is.na(ageInt), case_mask[row_number()]) %>% mutate(outcome = 1L)
  ctrls <- base %>% filter(!is.na(ageInt), status == 0, BBD_history == 0) %>% mutate(outcome = 0L)
  d <- bind_rows(cases, ctrls) %>% mutate(study = factor(study))
  round(exp(coef(glm(outcome ~ c + ageInt + study, d, family = binomial()))["c"]), 2)
}

mask_A <- dat$status == 0 & dat$BBD_history == 1 & !is.na(dat$BBD_history)
mask_B <- dat$BBD_history == 1 & !is.na(dat$BBD_history)
contrast <- map_dfr(c("CHEK2","ATM","BRCA2"), function(g) {
  tibble(Gene = g,
         `(A) cancer-free BBD [primary]` = or_def(mask_A, g),
         `(B) any-status BBD [conflated]` = or_def(mask_B, g))
})

cat("=== (1) Carrier depletion of the cancer-free BBD case group ===\n")
print(as.data.frame(depletion))
cat("\n=== (2) BBD OR under two case definitions (A primary vs B includes cancer) ===\n")
print(as.data.frame(contrast))
cat("\nMethod (C) — adjusting for cancer status — is not estimable: every case in\n")
cat("definition (A) is status==0, so cancer status has no variation. A collider\n")
cat("cannot be removed by adjustment. The cancer-free BBD ORs (<1 for CHEK2/ATM)\n")
cat("are selection artifacts, not protective effects.\n")

write_xlsx(list(carrier_depletion = depletion, case_definition_contrast = contrast),
           "outputs/tables/bbd_collider_check.xlsx")
message("\nSaved: outputs/tables/bbd_collider_check.xlsx")
