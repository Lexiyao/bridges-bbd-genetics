# =============================================================================
# Sensitivity Analysis — Familial vs Population-Based Ascertainment
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil Population Health Sciences — 2025/26
#
# PURPOSE
#   The BRIDGES consortium pools population-based and familial (high-risk
#   clinic) studies. Familial-ascertainment studies (e.g. HEBON, GC-HBOC,
#   KCONFAB) preferentially recruit pathogenic variant carriers, which can
#   inflate case-control odds ratios for moderate-penetrance genes.
#
#   This script quantifies that influence by re-estimating each truncating
#   variant association under three specifications:
#     (A) All studies + study adjustment   [= primary analysis]
#     (B) Population-based studies only     [familial studies removed]
#     (C) All studies + adjustment for family history (famHist)
#
#   A study is flagged "familial-ascertainment" if >50% of its participants
#   report a positive family history (famHist == 1). The two largest such
#   contributors (HEBON, GC-HBOC) are hereditary breast/ovarian cancer
#   cohorts by design, so the flag is robust to the exact threshold.
#
# INTERPRETATION
#   - Genes whose OR is stable across (A)-(C) are robust to ascertainment.
#   - Genes whose OR attenuates under (B) are partly driven by familial
#     enrichment and should be reported with that caveat.
#
# DATA GOVERNANCE
#   Data accessed under BCAC/BRIDGES approved secondary use agreement.
#   All outputs are non-disclosive summary statistics.
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c("tidyverse", "logistf", "writexl")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(tidyverse)
library(logistf)
library(writexl)

# ── 1. PATHS AND CONSTANTS ────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
MISS_FILE  <- "concept_807_zhang_bridges_missense.csv"

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

GENES <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
           "BARD1", "RAD51C", "RAD51D", "TP53")

FAMILIAL_FH_THRESHOLD <- 0.50  # study flagged familial if >50% report family history

# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
message("Loading phenotype file...")
pheno <- read_delim(PHENO_FILE, delim = "\t", show_col_types = FALSE,
                    na = c("", "NA", "888", "777", "999"))

message("Loading truncating variants file...")
trunc <- read_csv(TRUNC_FILE, show_col_types = FALSE, na = c("", "NA")) %>%
  select(-any_of("...1"))

message("Loading missense variants file...")
miss <- read_csv(MISS_FILE, show_col_types = FALSE, na = c("", "NA")) %>%
  select(-any_of("...1"))

# Build a BRIDGES_ID + 9-gene 0/1 carrier table for one variant source.
#   Truncating carrier  : GENE_truncating == 1
#   Missense carrier     : GENE_CADD.phred.01 == 1 (CADD phred >= 20, top 1%)
# Columns are suffixed (_t / _m) so both can coexist on one merged data frame.
build_carriers <- function(raw, pattern, suffix) {
  cols <- map_chr(GENES, function(g) {
    hit <- names(raw)[grepl(paste0("^", g, pattern), names(raw), ignore.case = TRUE)]
    if (length(hit) > 0) hit[1] else NA_character_
  })
  keep <- !is.na(cols)
  raw %>%
    select(BRIDGES_ID, all_of(cols[keep])) %>%
    rename(!!!setNames(cols[keep], paste0(GENES[keep], suffix))) %>%
    mutate(across(-BRIDGES_ID, ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))
}

carriers_t <- build_carriers(trunc, "_truncating$",   "_t")
carriers_m <- build_carriers(miss,  "_CADD\\.phred\\.01$", "_m")

dat <- pheno %>%
  filter(ethnicityClass == 1, !is.na(ageInt)) %>%
  left_join(carriers_t, by = "BRIDGES_ID") %>%
  left_join(carriers_m, by = "BRIDGES_ID") %>%
  mutate(across(c(ends_with("_t"), ends_with("_m")), ~ replace_na(.x, 0L)))

genes_found <- GENES  # all 9 present in both files (verified by data_QC.R)

# ── 3. FLAG FAMILIAL-ASCERTAINMENT STUDIES ────────────────────────────────────
fam_studies <- dat %>%
  group_by(study) %>%
  summarise(pct_fh = mean(famHist == 1, na.rm = TRUE), .groups = "drop") %>%
  filter(pct_fh > FAMILIAL_FH_THRESHOLD) %>%
  pull(study)

message("\nFamilial-ascertainment studies (>", 100 * FAMILIAL_FH_THRESHOLD,
        "% family history):")
message("  ", paste(fam_studies, collapse = ", "))

# ── 4. COHORTS ────────────────────────────────────────────────────────────────
controls <- dat %>% filter(status == 0)
dcis     <- dat %>% filter(status == 2, MorphologygroupIndex_corr == "Ductal")

message(sprintf("\nDCIS cases: %d (%d from familial studies, %d population)",
                nrow(dcis),
                sum(dcis$study %in% fam_studies),
                sum(!dcis$study %in% fam_studies)))

# ── 5. FIRTH ESTIMATOR UNDER A GIVEN SPECIFICATION ────────────────────────────
# vtype = "truncating" -> column GENE_t ; "missense" -> column GENE_m
fit_spec <- function(gene, spec, vtype) {
  cs <- dcis; ct <- controls; adjust_fh <- FALSE
  carrier_col <- paste0(gene, if (vtype == "missense") "_m" else "_t")

  if (spec == "population") {
    cs <- cs %>% filter(!study %in% fam_studies)
    ct <- ct %>% filter(!study %in% fam_studies)
  } else if (spec == "famHist_adjusted") {
    adjust_fh <- TRUE
  }

  d <- bind_rows(ct %>% mutate(outcome = 0L),
                 cs %>% mutate(outcome = 1L)) %>%
    mutate(carrier = .data[[carrier_col]], study = factor(study))

  if (adjust_fh) d <- d %>% filter(!is.na(famHist))

  n_cc <- sum(d$carrier == 1 & d$outcome == 1)
  n_ct <- sum(d$carrier == 1 & d$outcome == 0)

  # Skip if no carriers in cases — OR would be a prior-driven artifact
  if (n_cc == 0) {
    return(tibble(gene, vtype, spec, n_case = sum(d$outcome == 1),
                  carr_case = n_cc, carr_ctrl = n_ct,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, note = "0 carriers in cases"))
  }

  form <- if (adjust_fh) outcome ~ carrier + ageInt + study + famHist
          else            outcome ~ carrier + ageInt + study
  fit <- logistf(form, data = d, firth = TRUE, pl = TRUE, plconf = 2)
  i   <- which(names(coef(fit)) == "carrier")

  tibble(gene, vtype, spec, n_case = sum(d$outcome == 1),
         carr_case = n_cc, carr_ctrl = n_ct,
         OR      = exp(coef(fit)[i]),
         CI_low  = exp(fit$ci.lower[i]),
         CI_high = exp(fit$ci.upper[i]),
         p_value = fit$prob[i],
         note    = "")
}

# ── 6. RUN ALL GENES × 2 VARIANT TYPES × 3 SPECIFICATIONS ─────────────────────
message("\nRunning DCIS sensitivity (2 variant types × 3 specifications × ",
        length(genes_found), " genes)...")

specs  <- c("all_studies", "population", "famHist_adjusted")
vtypes <- c("truncating", "missense")
results <- map_dfr(vtypes, function(vt) {
  map_dfr(genes_found, function(g) {
    map_dfr(specs, function(s) {
      message(sprintf("  %-10s %-8s / %s", vt, g, s))
      fit_spec(g, s, vt)
    })
  })
})

# Apply BH-FDR within each variant type × specification (9 genes)
results <- results %>%
  group_by(vtype, spec) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# ── 7. WIDE COMPARISON TABLE ──────────────────────────────────────────────────
fmt_or <- function(or, lo, hi) {
  if_else(is.na(or), "—", sprintf("%.2f (%.2f–%.2f)", or, lo, hi))
}

wide <- results %>%
  mutate(cell = fmt_or(OR, CI_low, CI_high)) %>%
  select(gene, vtype, spec, cell) %>%
  pivot_wider(names_from = spec, values_from = cell) %>%
  left_join(
    results %>% filter(spec == "all_studies") %>%
      select(gene, vtype, carr_case_all = carr_case, carr_ctrl_all = carr_ctrl),
    by = c("gene", "vtype")
  ) %>%
  mutate(gene  = factor(gene, levels = GENES),
         vtype = factor(vtype, levels = c("truncating", "missense"))) %>%
  arrange(vtype, gene) %>%
  select(Variant = vtype,
         Gene = gene,
         `Carriers (case)`  = carr_case_all,
         `Carriers (ctrl)`  = carr_ctrl_all,
         `(A) All studies [PRIMARY]` = all_studies,
         `(B) Population only`       = population,
         `(C) famHist-adjusted`      = famHist_adjusted)

cat("\n", strrep("=", 78), "\n",
    "DCIS — SENSITIVITY TO FAMILIAL ASCERTAINMENT (truncating + missense)\n",
    strrep("=", 78), "\n", sep = "")
print(as.data.frame(wide))

# ── 8. EXPORT ─────────────────────────────────────────────────────────────────
out_path <- "outputs/tables/sensitivity_familial_population.xlsx"
write_xlsx(
  list(
    "Comparison"  = wide,
    "Raw_long"    = results,
    "Fam_studies" = tibble(familial_study = fam_studies)
  ),
  path = out_path
)
message(sprintf("\nExported: %s", out_path))
message("Done.")
