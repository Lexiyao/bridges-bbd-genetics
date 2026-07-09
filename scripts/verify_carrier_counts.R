# =============================================================================
# Independent Carrier-Count Verification (audit / provenance record)
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — University of Cambridge MPhil
#
# PURPOSE
#   Re-count every "Carriers (cases)" / "Carriers (controls)" value reported in
#   the result tables, using a code path that is INDEPENDENT of the analysis
#   scripts: cohort IDs are derived directly from the phenotype file, and
#   carriers are tallied straight from the raw genotype files for those IDs
#   (dominant model: value >= 1, which also captures multi-allelic value = 2).
#   The recounts are then compared against the published result tables. Any
#   discrepancy is flagged. This provides an auditable record that the carrier
#   counts — including the small counts for rare genes (TP53, RAD51D, BARD1) —
#   are correct and reflect true rarity, not missing data or a coding error.
#
# OUTPUT
#   outputs/QC/carrier_count_verification.xlsx   (full recount-vs-table table)
#   outputs/QC/carrier_count_verification_<ts>.txt (plain-text log with verdict)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(writexl)
})

dir.create("outputs/QC", recursive = TRUE, showWarnings = FALSE)
ts  <- format(Sys.time(), "%Y%m%d_%H%M%S")
log <- file(sprintf("outputs/QC/carrier_count_verification_%s.txt", ts), open = "wt")
say <- function(...) { m <- paste0(...); message(m); writeLines(m, log) }

GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")

# ── Load raw data ─────────────────────────────────────────────────────────────
pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                    show_col_types = FALSE, na = c("", "NA", "888", "777", "999"))
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv",
                  show_col_types = FALSE, na = c("", "NA"))
miss  <- read_csv("concept_807_zhang_bridges_missense.csv",
                  show_col_types = FALSE, na = c("", "NA"))

# Independent carrier count: dominant model on the raw genotype value for a set of IDs
count_carriers <- function(ids, gtab, gene, pattern) {
  col <- names(gtab)[grepl(paste0("^", gene, pattern, "$"), names(gtab))][1]
  v <- gtab[[col]][gtab$BRIDGES_ID %in% ids]
  sum(!is.na(v) & as.numeric(v) >= 1)
}

# ── Cohort ID sets (derived directly from phenotype) ──────────────────────────
ids <- list(
  ctrl = pheno %>% filter(status==0, ethnicityClass==1, !is.na(ageInt)) %>% pull(BRIDGES_ID),
  bbdc = pheno %>% filter(status==0, ethnicityClass==1, !is.na(ageInt), BBD_history==1) %>% pull(BRIDGES_ID),
  bbdk = pheno %>% filter(status==0, ethnicityClass==1, !is.na(ageInt), BBD_history==0) %>% pull(BRIDGES_ID),
  dcis = pheno %>% filter(status==2, ethnicityClass==1, !is.na(ageInt), MorphologygroupIndex_corr=="Ductal")  %>% pull(BRIDGES_ID),
  lcis = pheno %>% filter(status==2, ethnicityClass==1, !is.na(ageInt), MorphologygroupIndex_corr=="Lobular") %>% pull(BRIDGES_ID)
)

# ── Published tables ──────────────────────────────────────────────────────────
rd  <- function(f,s) as.data.frame(read_excel(f,s))
tbl <- list(
  bbd_trunc  = rd("outputs/tables/BBD_truncating_FINAL.xlsx",         "BBD vs Controls"),
  dcis_trunc = rd("outputs/tables/DCIS_LCIS_truncating_results.xlsx", "DCIS_truncating"),
  lcis_trunc = rd("outputs/tables/DCIS_LCIS_truncating_results.xlsx", "LCIS_truncating"),
  dcis_miss  = rd("outputs/tables/missense_results_FINAL.xlsx",       "DCIS_missense")
)

# (panel, variant pattern, case-id, ctrl-id, table)
panels <- list(
  list("BBD truncating",  "_truncating",        "bbdc","bbdk","bbd_trunc"),
  list("DCIS truncating", "_truncating",        "dcis","ctrl","dcis_trunc"),
  list("LCIS truncating", "_truncating",        "lcis","ctrl","lcis_trunc"),
  list("DCIS missense",   "_CADD\\.phred\\.01", "dcis","ctrl","dcis_miss")
)

rows <- list()
for (p in panels) {
  panel <- p[[1]]; pat <- p[[2]]; gt <- if (grepl("missense", panel)) miss else trunc
  tt <- tbl[[p[[5]]]]
  for (g in GENES) {
    rc <- count_carriers(ids[[p[[3]]]], gt, g, pat)
    rk <- count_carriers(ids[[p[[4]]]], gt, g, pat)
    tr <- tt[tt$Gene == g, ]
    tc <- as.integer(tr[["Carriers (cases)"]]); tk <- as.integer(tr[["Carriers (controls)"]])
    rows[[length(rows)+1]] <- tibble(
      panel = panel, gene = g,
      recount_cases = rc, table_cases = tc, cases_match = rc == tc,
      recount_ctrls = rk, table_ctrls = tk, ctrls_match = rk == tk)
  }
}
out <- bind_rows(rows)
n_mismatch <- sum(!out$cases_match | !out$ctrls_match)

say("Independent Carrier-Count Verification")
say("Run: ", format(Sys.time()))
say(strrep("=", 60))
say("Method: cohort IDs from phenotype; carriers tallied directly from raw")
say("genotype files (dominant, value>=1). Independent of analysis scripts.")
say("")
for (p in unique(out$panel)) {
  say("-- ", p, " --")
  sub <- out %>% filter(panel == p)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i,]
    say(sprintf("  %-7s cases %d/%d %s   ctrls %d/%d %s", r$gene,
        r$recount_cases, r$table_cases, ifelse(r$cases_match,"OK","MISMATCH"),
        r$recount_ctrls, r$table_ctrls, ifelse(r$ctrls_match,"OK","MISMATCH")))
  }
}
say("")
say(strrep("=", 60))
say(sprintf("VERDICT: %d of %d carrier counts re-verified; %d mismatch(es).",
            2*nrow(out), 2*nrow(out), n_mismatch))
say(if (n_mismatch == 0)
      "ALL carrier counts independently confirmed correct (incl. small counts)."
    else "*** MISMATCH FOUND — investigate before reporting. ***")
say("")
say("Note: the truncating/missense genotype files contain no missing calls")
say("(verified in data_QC.R), so small counts reflect true rarity of these")
say("genes (e.g. TP53, RAD51D), not missing data.")

write_xlsx(list(carrier_count_check = out),
           "outputs/QC/carrier_count_verification.xlsx")
say("\nSaved: outputs/QC/carrier_count_verification.xlsx")
close(log)
message("Done.")
