# =============================================================================
# KARMA-updated BBD status — sensitivity re-run for the BBD genetic analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE
#   After the primary analysis was finalised, the KARMA study supplied an
#   updated version of its BBD-status data (the original KARMA BBD field was
#   near-empty, contributing ~1 case to the analysis cohort). This script
#   quantifies the robustness of the BBD truncating-variant results to that
#   update, comparing three cohorts (all filtered to the primary analysis
#   frame: status==0 & ethnicityClass==1 & non-missing BBD_history & ageInt):
#     (PRIMARY)  pheno as-is, KARMA old coding      ~1765 cases / 12715 controls
#     (EXCLUDED) PRIMARY minus ALL study=="KARMA"   ~1764 cases /  7238 controls
#     (UPDATED)  KARMA BBD status reclassified      ~2045 cases / 12400 controls
#   The dissertation reports EXCLUDED (drop-all-KARMA) vs UPDATED as the two
#   robustness definitions; PRIMARY is the main analysis. Note EXCLUDED is NOT
#   PRIMARY: the two differ by the ~5,478 KARMA participants (BRCA2 overall CI
#   upper 6.64 vs 6.63 respectively — distinct cohorts, not a rounding artefact).
#
# DATA AVAILABILITY
#   Individual-level BRIDGES/BCAC data and the KARMA BBD-status update are
#   accessed under an approved data agreement and are NOT distributed with this
#   repository. The input filenames below are placeholders; the script is
#   provided for methodological transparency and is run from the repository root.
#   Aggregate cells with <5 carriers are non-disclosive and are not released.
#
# INPUTS (not distributed)
#   concept_807_zhang_bridges_pheno_v17.txt    phenotype table
#   concept_807_zhang_bridges_truncating.csv   per-gene truncating carrier flags
#   karma_bbd_update.xlsx    KARMA BBD-status update (one ';'-packed column)
#   karma_id_map.txt         uniqueid -> BCAC_ID crosswalk
#
# OUTPUT
#   outputs/tables/karma_bbd_sensitivity.xlsx
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(readxl); library(logistf); library(writexl)})
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

GENES  <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")
PANELS <- c("Overall BBD", "Non-proliferative", "Proliferative")
NA_CODES <- c("", "NA", "888", "777", "999")

KARMA_XLSX  <- "karma_bbd_update.xlsx"   # KARMA BBD-status update; not distributed
KARMA_MATCH <- "karma_id_map.txt"        # uniqueid -> BCAC_ID crosswalk; not distributed

pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                    show_col_types = FALSE, na = NA_CODES)
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv", show_col_types = FALSE, na = c("", "NA"))
trg   <- trunc %>% select(BRIDGES_ID, all_of(paste0(GENES, "_truncating")))
names(trg) <- c("BRIDGES_ID", GENES)

# ---- build the KARMA-UPDATED pheno (overwrite 6 BBD fields, KARMA rows only) ----
raw  <- read_excel(KARMA_XLSX, col_names = FALSE)
hdr  <- strsplit(raw[[1]][1], ";")[[1]]
recs <- lapply(raw[[1]][-1], function(x) { r <- strsplit(x, ";")[[1]]; length(r) <- length(hdr); r })
kd   <- as.data.frame(do.call(rbind, recs), stringsAsFactors = FALSE); names(kd) <- hdr
bbd6 <- c("BBD_history","BBD_number","BBD_type1","BBD_type2","BBD_type3","BBD_type4")
kd[bbd6] <- lapply(kd[bbd6], function(v) ifelse(v %in% NA_CODES, NA, v))

match_tbl <- read_delim(KARMA_MATCH, delim = "\t", show_col_types = FALSE, na = character(0)) %>%
  select(uniqueid, BCAC_ID) %>% filter(!is.na(BCAC_ID) & BCAC_ID != "")
kd2 <- kd %>% inner_join(match_tbl, by = "uniqueid") %>% select(BCAC_ID, all_of(bbd6))
kd2 <- kd2[!duplicated(kd2$BCAC_ID), ]

phenoU <- pheno
hit <- phenoU$study == "KARMA" & phenoU$BCAC_ID %in% kd2$BCAC_ID
idx <- match(phenoU$BCAC_ID[hit], kd2$BCAC_ID)
for (f in bbd6) phenoU[[f]][hit] <- kd2[[f]][idx]

# ---- cohort builders ----
frame <- function(d) d %>%
  left_join(trg, by = "BRIDGES_ID") %>%
  mutate(ageInt = suppressWarnings(as.numeric(ageInt))) %>%
  filter(status == "0", ethnicityClass == "1", !is.na(ageInt))

cohorts <- list(
  PRIMARY  = frame(pheno),
  excluded = frame(pheno) %>% filter(study != "KARMA"),
  updated  = frame(phenoU)
)

# ---- one Firth fit ----
fit_cell <- function(d, gene, panel) {
  d$carrier <- as.integer(replace_na(as.numeric(d[[gene]]), 0) >= 1)
  d$y <- switch(panel,
    "Overall BBD"       = ifelse(d$BBD_history == "1", 1L, ifelse(d$BBD_history == "0", 0L, NA)),
    "Non-proliferative" = ifelse(d$BBD_history == "1" & d$BBD_type1 == "1", 1L, ifelse(d$BBD_history == "0", 0L, NA)),
    "Proliferative"     = ifelse(d$BBD_history == "1" & d$BBD_type1 == "2", 1L, ifelse(d$BBD_history == "0", 0L, NA)))
  d <- d %>% filter(!is.na(y))
  out <- tibble(gene = gene, panel = panel,
                carrier_cases = sum(d$carrier == 1 & d$y == 1),
                carrier_controls = sum(d$carrier == 1 & d$y == 0),
                cases = sum(d$y == 1), controls = sum(d$y == 0),
                OR = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_)
  if (out$carrier_cases > 0) {
    f <- try(logistf(y ~ carrier + ageInt + factor(study), data = d, pl = TRUE), silent = TRUE)
    if (!inherits(f, "try-error")) {
      ci <- confint(f)["carrier", ]
      out$OR <- exp(coef(f)["carrier"]); out$ci_lo <- exp(ci[1]); out$ci_hi <- exp(ci[2]); out$p <- f$prob["carrier"]
    }
  }
  out
}

run_arm <- function(d, tag) {
  res <- map_dfr(PANELS, function(pn) map_dfr(GENES, ~ fit_cell(d, .x, pn)))
  res %>% group_by(panel) %>% mutate(q_fdr = p.adjust(p, "BH")) %>% ungroup() %>% mutate(analysis = tag)
}

results <- bind_rows(run_arm(cohorts$PRIMARY, "primary"),
                     run_arm(cohorts$excluded, "excluded"),
                     run_arm(cohorts$updated, "updated"))

write_xlsx(list("KARMA_sensitivity" = results), "outputs/tables/karma_bbd_sensitivity.xlsx")

# ---- console summary (headline cells) ----
show <- results %>%
  filter(panel == "Overall BBD", gene %in% c("BRCA2", "BARD1", "CHEK2")) %>%
  mutate(across(c(OR, ci_lo, ci_hi), ~ round(., 2)), p = round(p, 3), q_fdr = round(q_fdr, 3)) %>%
  select(analysis, gene, carrier_cases, cases, controls, OR, ci_lo, ci_hi, p, q_fdr) %>%
  arrange(gene, analysis)
cat("\n=== Overall BBD, headline genes (primary / excluded / updated) ===\n")
print(as.data.frame(show))
cat("\nWrote outputs/tables/karma_bbd_sensitivity.xlsx\n")
