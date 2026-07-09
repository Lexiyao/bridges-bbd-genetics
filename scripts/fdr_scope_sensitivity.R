# =============================================================================
# FDR Correction Scope — Sensitivity to the multiple-testing family definition
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE (addresses supervisor point #5: "FDR correction scope")
#   The primary analysis applies Benjamini-Hochberg WITHIN each outcome ×
#   variant-class family (9 genes each: 6 families = BBD/DCIS/LCIS × trunc/miss).
#   A reviewer may argue the family should be broader (all 54 primary tests
#   jointly) or that a more stringent control is warranted. This script
#   recomputes significance under three scopes and shows which findings are
#   robust to the choice:
#     (1) Per-family BH (9 tests)      = PRIMARY, as reported
#     (2) Global BH across 54 tests    = single joint family
#     (3) Global Bonferroni (54 tests) = conservative FWER reference
#
#   Subtype analyses (non-proliferative / proliferative BBD) are EXPLORATORY
#   and excluded from the primary multiple-testing family (documented as such).
#
# Reads only existing result tables — no model refitting.
# DATA GOVERNANCE: BCAC/BRIDGES approved secondary use; non-disclosive output.
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
})

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ── 1. PULL FULL-PRECISION p-VALUES FROM THE THREE RESULT FILES ────────────────
raw <- bind_rows(
  read_excel("outputs/tables/BBD_truncating_FINAL.xlsx",        "Raw results"),
  read_excel("outputs/tables/DCIS_LCIS_truncating_results.xlsx","Raw_results"),
  read_excel("outputs/tables/missense_results_FINAL.xlsx",      "All_raw")
) %>%
  select(analysis, gene, OR, p_value)

# ── 2. TAG OUTCOME, VARIANT CLASS, AND PRIMARY-vs-SUBTYPE ──────────────────────
PRIMARY_ANALYSES <- c(
  "BBD vs Controls", "DCIS vs Controls (truncating)", "LCIS vs Controls (truncating)",
  "BBD vs Controls (missense)", "DCIS vs Controls (missense)", "LCIS vs Controls (missense)"
)

tagged <- raw %>%
  mutate(
    variant = if_else(str_detect(analysis, "missense"), "missense", "truncating"),
    outcome = case_when(
      str_detect(analysis, "^DCIS")            ~ "DCIS",
      str_detect(analysis, "^LCIS")            ~ "LCIS",
      str_detect(analysis, "Non-prolif")       ~ "BBD non-prolif",
      str_detect(analysis, "Prolif")           ~ "BBD prolif",
      str_detect(analysis, "^BBD vs Controls") ~ "BBD",
      TRUE                                     ~ analysis),
    family  = paste(outcome, variant),
    is_primary = analysis %in% PRIMARY_ANALYSES
  ) %>%
  filter(!is.na(p_value))

primary <- tagged %>% filter(is_primary)
message(sprintf("Primary tests: %d (expect 54 = 6 families × 9 genes)", nrow(primary)))

# ── 3. THREE CORRECTION SCOPES ────────────────────────────────────────────────
scored <- primary %>%
  group_by(family) %>%
  mutate(q_perfamily_BH = p.adjust(p_value, method = "BH")) %>%   # (1) PRIMARY
  ungroup() %>%
  mutate(
    q_global_BH   = p.adjust(p_value, method = "BH"),             # (2) 54 jointly
    p_bonferroni  = p.adjust(p_value, method = "bonferroni"),     # (3) FWER ref
    sig_perfamily = q_perfamily_BH < 0.05,
    sig_globalBH  = q_global_BH   < 0.05,
    sig_bonf      = p_bonferroni  < 0.05
  )

# ── 4. REPORT: every finding significant under ANY scope ──────────────────────
report <- scored %>%
  filter(sig_perfamily | sig_globalBH | sig_bonf) %>%
  arrange(p_value) %>%
  transmute(
    outcome, variant, gene,
    OR = sprintf("%.2f", OR),
    p_value = if_else(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    `q per-family BH (PRIMARY)` = sprintf("%.3f", q_perfamily_BH),
    `q global BH (54)`          = sprintf("%.3f", q_global_BH),
    `p Bonferroni (54)`         = if_else(p_bonferroni < 0.001, "<0.001", sprintf("%.3f", p_bonferroni)),
    survives = case_when(
      sig_bonf     ~ "ALL scopes incl. Bonferroni — strongest",
      sig_globalBH ~ "per-family + global BH (not Bonferroni)",
      TRUE         ~ "per-family BH only — scope-sensitive"
    ))

cat("\n", strrep("=", 96), "\n",
    "FDR SCOPE SENSITIVITY — findings significant under at least one correction scope\n",
    "(", nrow(primary), " primary tests: 6 outcome×variant families × 9 genes)\n",
    strrep("=", 96), "\n", sep = "")
as.data.frame(report) %>% print(right = FALSE)

cat("\nCounts of FDR-significant primary findings by scope:\n")
cat(sprintf("  (1) per-family BH (PRIMARY): %d\n", sum(scored$sig_perfamily)))
cat(sprintf("  (2) global BH (54 jointly) : %d\n", sum(scored$sig_globalBH)))
cat(sprintf("  (3) Bonferroni (54)        : %d\n", sum(scored$sig_bonf)))

# ── 5. EXPORT ─────────────────────────────────────────────────────────────────
write_xlsx(
  list("Scope_report" = as.data.frame(report),
       "All_primary_scored" = scored %>% as.data.frame()),
  "outputs/tables/fdr_scope_sensitivity.xlsx")
message("\nExported: outputs/tables/fdr_scope_sensitivity.xlsx")
