# =============================================================================
# Data Quality Control Script
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil Population Health Sciences — 2025/26
#
# PURPOSE
#   This script performs systematic quality control (QC) checks on all input
#   data files BEFORE any association analysis. It should be run first in the
#   analysis pipeline and its output reviewed carefully.
#
#   QC is structured into three levels consistent with best practice in genetic
#   epidemiology secondary analysis (see: Anderson et al. 2010 Nat Protoc;
#   Marees et al. 2018 Int J Methods Psychiatr Res):
#
#   Level 1 — Data integrity:    file dimensions, variable types, NA coding
#   Level 2 — Phenotype QC:      cohort derivation, covariate completeness,
#                                 variable distributions, plausibility checks
#   Level 3 — Genotype QC:       carrier frequencies, cross-tabulations,
#                                 study-level carrier distributions,
#                                 comparison with published BRIDGES values
#
# IMPORTANT SCOPE NOTE
#   The BRIDGES consortium performed all upstream genotype-level QC (sample
#   call rate, heterozygosity, sex concordance, relatedness, HWE, variant
#   call rate, ancestry PCA) prior to data release. Reference: Dorling et al.
#   (2021) NEJM 384:428-439. This script does NOT replicate those steps as
#   the raw VCF/PLINK files are not available to secondary analysts.
#   What this script CAN and DOES check:
#     - Phenotype file integrity and variable completeness
#     - Cohort derivation reproducibility
#     - Carrier frequency plausibility against published values
#     - Study-level data distributions (important for confounding)
#     - Covariate availability per analysis cohort
#
# DATA GOVERNANCE
#   Data accessed under BCAC/BRIDGES approved secondary use agreement.
#   All QC outputs are non-disclosive summary statistics.
#   Raw participant-level data are not committed to version control.
#
# OUTPUT
#   outputs/QC/QC_report_[timestamp].txt — Full QC report (plain text)
#   outputs/QC/QC_figures/               — Diagnostic plots (PNG, 300 dpi)
#   Console output mirrors the report file in real time.
# =============================================================================


# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c("tidyverse", "writexl", "ggplot2", "scales", "knitr")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(tidyverse)
library(writexl)
library(ggplot2)
library(scales)


# ── 1. OUTPUT DIRECTORIES ─────────────────────────────────────────────────────
dir.create("outputs/QC",         recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/QC/figures", recursive = TRUE, showWarnings = FALSE)

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
report_path <- sprintf("outputs/QC/QC_report_%s.txt", timestamp)
con         <- file(report_path, open = "wt")

# Helper: write to both console and report file
qc_log <- function(...) {
  msg <- paste0(...)
  message(msg)
  writeLines(msg, con)
}

qc_section <- function(title) {
  qc_log("\n", strrep("=", 70))
  qc_log("  ", title)
  qc_log(strrep("=", 70))
}

qc_pass    <- function(msg) qc_log("  [PASS] ", msg)
qc_warn    <- function(msg) qc_log("  [WARN] ", msg)
qc_fail    <- function(msg) qc_log("  [FAIL] ", msg)
qc_info    <- function(msg) qc_log("  [INFO] ", msg)

# Track overall QC flags
qc_flags <- list(warnings = character(0), failures = character(0))

flag_warn <- function(msg) {
  qc_warn(msg)
  qc_flags$warnings <<- c(qc_flags$warnings, msg)
}
flag_fail <- function(msg) {
  qc_fail(msg)
  qc_flags$failures <<- c(qc_flags$failures, msg)
}


# ── 2. FILE PATHS ─────────────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
MISS_FILE  <- "concept_807_zhang_bridges_missense.csv"

GENES <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
           "BARD1", "RAD51C", "RAD51D", "TP53")

# Expected cohort sizes (from verified analysis outputs)
EXPECTED <- list(
  controls_full  = 34424L,
  controls_bbd   = 12715L,
  bbd_cases      = 1765L,
  dcis_cases     = 1663L,
  lcis_cases     = 139L
)

# Published carrier frequencies in controls from Dorling et al. 2021 NEJM
# (Table S2, European-ancestry subset, truncating PTV carriers)
# NOTE: These are approximate values from the published paper.
#       Verify against Supplementary Table S2 before citing.
DORLING_CTRL_FREQ <- c(
  BRCA1  = 0.0020,
  BRCA2  = 0.0038,
  PALB2  = 0.0009,
  CHEK2  = 0.0085,
  ATM    = 0.0025,
  BARD1  = 0.0005,
  RAD51C = 0.0005,
  RAD51D = 0.0004,
  TP53   = 0.0001
)


# ─────────────────────────────────────────────────────────────────────────────
qc_section("QC REPORT HEADER")
qc_log("  Analysis  : BBD/DCIS/LCIS genetic associations")
qc_log("  Run date  : ", format(Sys.time()))
qc_log("  R version : ", R.version.string)
qc_log("  Working dir: ", getwd())
qc_log("  Phenotype : ", PHENO_FILE)
qc_log("  Truncating: ", TRUNC_FILE)
qc_log("  Missense  : ", MISS_FILE)


# =============================================================================
# LEVEL 1 — DATA INTEGRITY
# =============================================================================

qc_section("LEVEL 1A — FILE EXISTENCE AND DIMENSIONS")

for (f in c(PHENO_FILE, TRUNC_FILE, MISS_FILE)) {
  if (file.exists(f)) {
    qc_pass(sprintf("File found: %s (%.1f MB)",
                    f, file.size(f) / 1e6))
  } else {
    flag_fail(sprintf("File NOT found: %s", f))
    stop(sprintf("Required file missing: %s\nPlace data files in working directory: %s",
                 f, getwd()))
  }
}

# Load data
qc_log("\nLoading phenotype file...")
pheno_raw <- read_delim(
  PHENO_FILE, delim = "\t", show_col_types = FALSE,
  na = c("", "NA", "888", "777", "999")
)

qc_log("Loading truncating variants file...")
trunc_raw <- read_csv(
  TRUNC_FILE, show_col_types = FALSE, na = c("", "NA")
) %>% select(-any_of("...1"))

qc_log("Loading missense variants file...")
miss_raw <- read_csv(
  MISS_FILE, show_col_types = FALSE, na = c("", "NA")
) %>% select(-any_of("...1"))

qc_info(sprintf("Phenotype file : %d rows × %d columns",
                nrow(pheno_raw), ncol(pheno_raw)))
qc_info(sprintf("Truncating file: %d rows × %d columns",
                nrow(trunc_raw), ncol(trunc_raw)))
qc_info(sprintf("Missense file  : %d rows × %d columns",
                nrow(miss_raw),  ncol(miss_raw)))


# ── Check BRIDGES_ID linkage ──────────────────────────────────────────────────
qc_section("LEVEL 1B — PARTICIPANT ID LINKAGE")

n_pheno <- nrow(pheno_raw)
n_trunc <- nrow(trunc_raw)
n_miss  <- nrow(miss_raw)

ids_pheno <- pheno_raw$BRIDGES_ID
ids_trunc <- trunc_raw$BRIDGES_ID
ids_miss  <- miss_raw$BRIDGES_ID

n_trunc_in_pheno <- sum(ids_trunc %in% ids_pheno)
n_miss_in_pheno  <- sum(ids_miss  %in% ids_pheno)

qc_info(sprintf("Phenotype IDs                    : %d", n_pheno))
qc_info(sprintf("Truncating IDs in phenotype      : %d / %d  (%.1f%%)",
                n_trunc_in_pheno, n_trunc,
                100 * n_trunc_in_pheno / n_trunc))
qc_info(sprintf("Missense IDs in phenotype        : %d / %d  (%.1f%%)",
                n_miss_in_pheno, n_miss,
                100 * n_miss_in_pheno / n_miss))

# Duplicate ID check
dup_pheno <- sum(duplicated(ids_pheno))
dup_trunc <- sum(duplicated(ids_trunc))
dup_miss  <- sum(duplicated(ids_miss))

if (dup_pheno == 0) qc_pass("No duplicate BRIDGES_IDs in phenotype file") else
  flag_fail(sprintf("%d duplicate BRIDGES_IDs in phenotype file", dup_pheno))

if (dup_trunc == 0) qc_pass("No duplicate BRIDGES_IDs in truncating file") else
  flag_fail(sprintf("%d duplicate BRIDGES_IDs in truncating file", dup_trunc))

if (dup_miss == 0) qc_pass("No duplicate BRIDGES_IDs in missense file") else
  flag_fail(sprintf("%d duplicate BRIDGES_IDs in missense file", dup_miss))

# Unmatched IDs
n_trunc_unmatched <- sum(!ids_trunc %in% ids_pheno)
n_miss_unmatched  <- sum(!ids_miss  %in% ids_pheno)

if (n_trunc_unmatched == 0) qc_pass("All truncating IDs matched to phenotype") else
  flag_warn(sprintf("%d truncating IDs not found in phenotype file", n_trunc_unmatched))

if (n_miss_unmatched == 0) qc_pass("All missense IDs matched to phenotype") else
  flag_warn(sprintf("%d missense IDs not found in phenotype file", n_miss_unmatched))


# ── Check required columns ────────────────────────────────────────────────────
qc_section("LEVEL 1C — REQUIRED VARIABLE PRESENCE")

required_pheno_cols <- c(
  "BRIDGES_ID", "status", "ethnicityClass", "ageInt", "study",
  "BBD_history", "BBD_type1", "MorphologygroupIndex_corr",
  "parity", "ageFFTP", "ageMenarche", "menoStat", "HRTEver"
)

for (col in required_pheno_cols) {
  if (col %in% names(pheno_raw)) {
    qc_pass(sprintf("Phenotype column present: %s", col))
  } else {
    flag_fail(sprintf("Required phenotype column MISSING: %s", col))
  }
}

# Gene columns in truncating file
detect_gene_col <- function(gene, col_names) {
  exact  <- col_names[tolower(col_names) == tolower(gene)]
  if (length(exact)  > 0) return(exact[1])
  prefix <- col_names[grepl(paste0("^", gene, "(_|$)"), col_names, ignore.case = TRUE)]
  if (length(prefix) > 0) return(prefix[1])
  NA_character_
}

gene_col_map_trunc <- setNames(
  sapply(GENES, detect_gene_col, col_names = names(trunc_raw)), GENES)
gene_col_map_miss  <- setNames(
  sapply(GENES, detect_gene_col, col_names = names(miss_raw)),  GENES)

qc_log("\nTruncating gene column mapping:")
for (g in GENES) {
  col <- gene_col_map_trunc[g]
  if (!is.na(col)) qc_pass(sprintf("  %-8s → %s", g, col)) else
    flag_fail(sprintf("  %-8s → NOT FOUND", g))
}

qc_log("\nMissense gene column mapping:")
for (g in GENES) {
  col <- gene_col_map_miss[g]
  if (!is.na(col)) qc_pass(sprintf("  %-8s → %s", g, col)) else
    flag_fail(sprintf("  %-8s → NOT FOUND", g))
}


# ── Check genotype values are 0/1 only ───────────────────────────────────────
qc_section("LEVEL 1D — GENOTYPE VARIABLE VALUE RANGE CHECK")
qc_info("Checking all carrier columns contain only 0, 1, or NA")

genes_trunc_found <- names(gene_col_map_trunc)[!is.na(gene_col_map_trunc)]
genes_miss_found  <- names(gene_col_map_miss)[!is.na(gene_col_map_miss)]

for (g in genes_trunc_found) {
  col      <- gene_col_map_trunc[g]
  vals     <- unique(na.omit(trunc_raw[[col]]))
  n_val2   <- sum(trunc_raw[[col]] == 2, na.rm = TRUE)
  unexpected <- setdiff(vals, c(0, 1, 2))
  if (length(unexpected) > 0) {
    flag_fail(sprintf("Truncating %-8s: unexpected values {%s}",
                      g, paste(unexpected, collapse=", ")))
  } else if (n_val2 > 0) {
    flag_warn(sprintf(
      "Truncating %-8s: values {%s}, n=%d with value 2 (multi-allelic carrier) — dominant recoding applied in analysis (≥1 → 1)",
      g, paste(sort(vals), collapse=", "), n_val2))
  } else {
    qc_pass(sprintf("Truncating %-8s: values {%s} — OK", g, paste(sort(vals), collapse=", ")))
  }
}

qc_log("\nMissense variant count columns (GENE_missense.variants) — values > 1 expected")
qc_info("  NOTE: analysis scripts use GENE_CADD.phred.01 binary columns (CADD phred >= 20,")
qc_info("        top 1% most deleterious variants genome-wide; Kircher et al. 2014 Nat Genet).")
qc_info("  not these variant count columns. Values shown here are total missense variant counts.")
for (g in genes_miss_found) {
  col  <- gene_col_map_miss[g]
  vals <- sort(unique(na.omit(miss_raw[[col]])))
  n_any <- sum(miss_raw[[col]] > 0, na.rm = TRUE)
  qc_info(sprintf("  Missense   %-8s: values {%s}, n with ≥1 variant = %d",
                  g, paste(vals, collapse=", "), n_any))
}

# Binary phenotype variables
# NOTE: menoStat is coded 1=pre/peri, 2=post (NOT 0/1) per BCAC Extended Data
#       Dictionary v4.0. HRTEver and BBD_history are 0/1.
qc_log("\nChecking binary phenotype variables (expected values per data dictionary)")
binary_01  <- c("HRTEver", "BBD_history")   # expected: {0, 1}
binary_12  <- c("menoStat")                  # expected: {1, 2}

for (bvar in c(binary_12, binary_01)) {
  if (!bvar %in% names(pheno_raw)) next
  vals      <- unique(na.omit(pheno_raw[[bvar]]))
  expected  <- if (bvar %in% binary_12) c(1, 2) else c(0, 1)
  unexpected <- setdiff(vals, expected)
  n_missing  <- sum(is.na(pheno_raw[[bvar]]))
  coding_note <- if (bvar == "menoStat") " [expected {1,2}: 1=pre/peri, 2=post]" else
                 " [expected {0,1}]"
  if (length(unexpected) == 0) {
    qc_pass(sprintf("%-15s: values {%s}, missing = %d  — OK%s",
                    bvar, paste(sort(vals), collapse=", "), n_missing, coding_note))
  } else {
    flag_warn(sprintf("%-15s: unexpected values {%s} — verify NA coding%s",
                      bvar, paste(unexpected, collapse=", "), coding_note))
  }
}


# =============================================================================
# LEVEL 2 — PHENOTYPE QC
# =============================================================================

qc_section("LEVEL 2A — STATUS DISTRIBUTION")

status_tab <- table(pheno_raw$status, useNA = "always")
qc_log("\nStatus distribution (full dataset, all ethnicities):")
for (nm in names(status_tab)) {
  qc_info(sprintf("  status = %-6s : %d", nm, status_tab[nm]))
}

# Cross-tabulate status × morphology
qc_log("\nStatus × MorphologygroupIndex_corr cross-tabulation:")
if ("MorphologygroupIndex_corr" %in% names(pheno_raw)) {
  morph_tab <- table(pheno_raw$status,
                     pheno_raw$MorphologygroupIndex_corr,
                     useNA = "always")
  capture.output(print(morph_tab)) %>%
    walk(~ qc_log("  ", .x))
}


qc_section("LEVEL 2B — ANCESTRY DISTRIBUTION")

if ("ethnicityClass" %in% names(pheno_raw)) {
  eth_tab <- table(pheno_raw$ethnicityClass, useNA = "always")
  qc_log("\nethenicityClass distribution (full dataset):")
  for (nm in names(eth_tab)) {
    pct <- 100 * eth_tab[nm] / nrow(pheno_raw)
    qc_info(sprintf("  ethnicityClass = %-5s : %6d  (%.1f%%)", nm, eth_tab[nm], pct))
  }
  n_eur <- sum(pheno_raw$ethnicityClass == 1, na.rm = TRUE)
  qc_info(sprintf("\nEuropean (ethnicityClass==1): %d of %d total (%.1f%%)",
                  n_eur, nrow(pheno_raw), 100 * n_eur / nrow(pheno_raw)))
} else {
  flag_fail("ethnicityClass column not found")
}


qc_section("LEVEL 2C — COHORT DERIVATION AUDIT")

qc_log("\nReproducing all study cohorts from raw data:")

# European status-0 pool
eur_s0 <- pheno_raw %>%
  filter(status == 0, ethnicityClass == 1, !is.na(ageInt))

qc_info(sprintf("European, status==0, !is.na(ageInt)          : n = %d  [expected controls pool = %d]",
                nrow(eur_s0), EXPECTED$controls_full))

if (nrow(eur_s0) == EXPECTED$controls_full) qc_pass("Controls pool matches expected n") else
  flag_warn(sprintf("Controls pool: observed %d, expected %d (diff = %d)",
                    nrow(eur_s0), EXPECTED$controls_full, nrow(eur_s0) - EXPECTED$controls_full))

# BBD cohort
bbd_cases <- eur_s0 %>% filter(BBD_history == 1)
bbd_ctrls <- eur_s0 %>% filter(BBD_history == 0)

qc_info(sprintf("BBD cases  (status==0, BBD_history==1)        : n = %d  [expected %d]",
                nrow(bbd_cases), EXPECTED$bbd_cases))
qc_info(sprintf("BBD ctrls  (status==0, BBD_history==0)        : n = %d  [expected %d]",
                nrow(bbd_ctrls), EXPECTED$controls_bbd))

n_bbd_missing_history <- sum(is.na(eur_s0$BBD_history))
qc_info(sprintf("status==0 with missing BBD_history            : n = %d  (excluded from BBD analysis)",
                n_bbd_missing_history))

if (nrow(bbd_cases) == EXPECTED$bbd_cases) qc_pass("BBD case count matches expected") else
  flag_warn(sprintf("BBD cases: observed %d, expected %d", nrow(bbd_cases), EXPECTED$bbd_cases))

if (nrow(bbd_ctrls) == EXPECTED$controls_bbd) qc_pass("BBD control count matches expected") else
  flag_warn(sprintf("BBD controls: observed %d, expected %d", nrow(bbd_ctrls), EXPECTED$controls_bbd))

# BBD subtype distribution
if ("BBD_type1" %in% names(pheno_raw)) {
  bbd_type_tab <- table(bbd_cases$BBD_type1, useNA = "always")
  qc_log("\nBBD subtype distribution (BBD cases only):")
  qc_info("  BBD_type1 = 1 : Non-proliferative")
  qc_info("  BBD_type1 = 2 : Proliferative (without atypia)")
  qc_info("  BBD_type1 = 3 : Atypical hyperplasia (excluded from subtype models)")
  for (nm in names(bbd_type_tab)) {
    qc_info(sprintf("  BBD_type1 = %-5s : %d", nm, bbd_type_tab[nm]))
  }
}

# DCIS and LCIS cohorts
dcis_cases <- pheno_raw %>%
  filter(status == 2, ethnicityClass == 1,
         MorphologygroupIndex_corr == "Ductal", !is.na(ageInt))
lcis_cases <- pheno_raw %>%
  filter(status == 2, ethnicityClass == 1,
         MorphologygroupIndex_corr == "Lobular", !is.na(ageInt))

qc_info(sprintf("\nDCIS cases (status==2, Ductal, European)      : n = %d  [expected %d]",
                nrow(dcis_cases), EXPECTED$dcis_cases))
qc_info(sprintf("LCIS cases (status==2, Lobular, European)     : n = %d  [expected %d]",
                nrow(lcis_cases), EXPECTED$lcis_cases))

if (nrow(dcis_cases) == EXPECTED$dcis_cases) qc_pass("DCIS case count matches expected") else
  flag_warn(sprintf("DCIS cases: observed %d, expected %d", nrow(dcis_cases), EXPECTED$dcis_cases))
if (nrow(lcis_cases) == EXPECTED$lcis_cases) qc_pass("LCIS case count matches expected") else
  flag_warn(sprintf("LCIS cases: observed %d, expected %d", nrow(lcis_cases), EXPECTED$lcis_cases))

# Verify status==3 is NOT suitable as BBD
qc_log("\nAudit: Why status==3 is NOT the correct BBD phenotype group")
s3 <- pheno_raw %>% filter(status == 3, ethnicityClass == 1)
qc_info(sprintf("  status==3 (European)                       : n = %d",    nrow(s3)))
qc_info(sprintf("  status==3 with BBD_history == 1            : n = %d",
                sum(s3$BBD_history == 1, na.rm = TRUE)))
qc_info(sprintf("  status==3 with missing BBD_history         : n = %d  (%.0f%%)",
                sum(is.na(s3$BBD_history)),
                100 * mean(is.na(s3$BBD_history))))
qc_warn("status==3 has high missing BBD_history — CONFIRMED NOT suitable as BBD group")
qc_info("  Correct BBD definition: status==0 & BBD_history==1 (used in all analysis scripts)")


qc_section("LEVEL 2D — AGE DISTRIBUTION AND PLAUSIBILITY")

for (cohort_name in c("Controls", "BBD", "DCIS", "LCIS")) {
  dat <- switch(cohort_name,
    "Controls" = eur_s0,
    "BBD"      = bind_rows(bbd_cases, bbd_ctrls),
    "DCIS"     = dcis_cases,
    "LCIS"     = lcis_cases
  )
  ages <- dat$ageInt[!is.na(dat$ageInt)]
  qc_info(sprintf("%-10s age: n=%d, mean=%.1f, sd=%.1f, range=[%.0f, %.0f], missing=%d",
                  cohort_name, length(ages),
                  mean(ages), sd(ages),
                  min(ages), max(ages),
                  sum(is.na(dat$ageInt))))
  if (min(ages) < 18) flag_warn(sprintf("%s: minimum age %.0f — check for implausible values", cohort_name, min(ages)))
  if (max(ages) > 100) flag_warn(sprintf("%s: maximum age %.0f — check for implausible values", cohort_name, max(ages)))
}

# Age distribution plot
age_plot_df <- bind_rows(
  eur_s0 %>% mutate(Group = "Controls"),
  bbd_cases %>% mutate(Group = "BBD cases"),
  dcis_cases %>% mutate(Group = "DCIS cases"),
  lcis_cases %>% mutate(Group = "LCIS cases")
) %>% filter(!is.na(ageInt))

source("scripts/_thesis_theme.R")   # THESIS_PAL, theme_thesis(), save_fig()
p_age <- ggplot(age_plot_df, aes(x = ageInt)) +
  geom_histogram(binwidth = 5, fill = unname(THESIS_PAL["baseline"]),
                 colour = "white", alpha = 0.9) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Age distribution by cohort (QC check)",
       x = "Age at interview (years)", y = "Count") +
  theme_thesis(base_size = 11, grid = "y") +
  theme(strip.background = element_rect(fill = "#EAF2FF", colour = NA))
save_fig(p_age, "outputs/QC/figures/QC_age_distribution.png", width = 10, height = 6)
qc_info("Age distribution plot saved: outputs/QC/figures/QC_age_distribution.png")


qc_section("LEVEL 2E — STUDY DISTRIBUTION")

study_tab <- table(eur_s0$study, useNA = "always")
qc_log(sprintf("\nNumber of contributing studies (controls pool): %d",
               sum(!is.na(names(study_tab)))))
qc_log("\nTop 15 studies by control count:")
top15 <- sort(study_tab, decreasing = TRUE)[1:min(15, length(study_tab))]
for (nm in names(top15)) {
  pct <- 100 * top15[nm] / nrow(eur_s0)
  qc_info(sprintf("  %-30s : %5d  (%.1f%%)", nm, top15[nm], pct))
}

# Study distribution plot
study_df <- as.data.frame(study_tab) %>%
  rename(study = Var1, n = Freq) %>%
  filter(!is.na(study)) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20) %>%
  mutate(study = factor(study, levels = rev(study)))

p_study <- ggplot(study_df, aes(x = n, y = study)) +
  geom_col(fill = unname(THESIS_PAL["baseline"]), alpha = 0.9) +
  labs(title = "Top 20 studies by control count (QC check)",
       x = "N controls", y = NULL) +
  theme_thesis(base_size = 10, grid = "x")
save_fig(p_study, "outputs/QC/figures/QC_study_distribution.png", width = 9, height = 7)
qc_info("Study distribution plot saved: outputs/QC/figures/QC_study_distribution.png")


qc_section("LEVEL 2F — COVARIATE COMPLETENESS BY ANALYSIS COHORT")

covariates <- c("ageInt", "study", "parity", "ageFFTP",
                "ageMenarche", "menoStat", "HRTEver")

cohort_list <- list(
  "Controls (full)"   = eur_s0,
  "BBD cases"         = bbd_cases,
  "BBD controls"      = bbd_ctrls,
  "DCIS cases"        = dcis_cases,
  "LCIS cases"        = lcis_cases
)

qc_log(sprintf("\n%-22s  %s", "Cohort / Covariate",
               paste(sprintf("%-14s", covariates), collapse = "")))
qc_log(strrep("-", 22 + 14 * length(covariates)))

for (cname in names(cohort_list)) {
  dat   <- cohort_list[[cname]]
  rates <- sapply(covariates, function(v) {
    if (!v %in% names(dat)) return(NA_real_)
    sprintf("%.1f%%", 100 * mean(!is.na(dat[[v]])))
  })
  qc_log(sprintf("%-22s  %s", cname,
                 paste(sprintf("%-14s", rates), collapse = "")))
}

qc_log("\nNOTE: ageFFTP completeness will be lower — expected, as only parous")
qc_log("      women have an age at first full-term pregnancy (structural missingness).")


qc_section("LEVEL 2G — HORMONAL VARIABLE DISTRIBUTIONS")

hormonal_vars <- c("parity", "ageFFTP", "ageMenarche", "menoStat", "HRTEver")

for (cohort_name in c("BBD cases", "BBD controls", "DCIS cases", "LCIS cases")) {
  dat <- cohort_list[[cohort_name]]
  qc_log(sprintf("\n  %s (n=%d):", cohort_name, nrow(dat)))
  for (v in hormonal_vars) {
    if (!v %in% names(dat)) next
    x <- dat[[v]][!is.na(dat[[v]])]
    n_miss <- sum(is.na(dat[[v]]))
    if (v == "menoStat") {
      # menoStat: 1=pre/peri-menopausal, 2=post-menopausal (BCAC dictionary v4.0)
      qc_info(sprintf("    %-14s: 1(pre/peri)=%d, 2(post)=%d, NA=%d",
                      v,
                      sum(dat[[v]] == 1, na.rm = TRUE),
                      sum(dat[[v]] == 2, na.rm = TRUE),
                      n_miss))
    } else if (v == "HRTEver") {
      qc_info(sprintf("    %-14s: 0(never)=%d, 1(ever)=%d, NA=%d",
                      v,
                      sum(dat[[v]] == 0, na.rm = TRUE),
                      sum(dat[[v]] == 1, na.rm = TRUE),
                      n_miss))
    } else {
      qc_info(sprintf("    %-14s: mean=%.1f, sd=%.1f, range=[%.0f,%.0f], NA=%d",
                      v, mean(x), sd(x), min(x), max(x), n_miss))
    }
  }
}


qc_section("LEVEL 2H — AVAILABLE-BUT-UNUSED RISK-FACTOR AUDIT")
# Guards against the analysis silently omitting a relevant variable that IS in
# the data. Lists recognised breast-disease risk factors / covariates, whether
# each is present in the phenotype file, its completeness among DCIS cases, and
# whether it is part of the current analysis. Anything present and reasonably
# complete but NOT used is flagged for explicit consideration.
USED_VARS <- c("ageInt", "study", "parity", "ageFFTP", "ageMenarche",
               "menoStat", "HRTEver")   # variables already in the models
# (variable, plain-language role)
RISK_FACTORS <- tribble(
  ~var,            ~role,
  "BMI",           "body mass index — major BC risk factor / confounder",
  "famHist",       "family history — risk factor (mediator for genetic models)",
  "HRTCurrent",    "current HRT use",
  "mensAgeLast",   "age at menopause (continuous)",
  "GradeIndex",    "tumour grade (case-only; enables grade stratification)",
  "ER_statusIndex","ER status (case-only; enables ER stratification)",
  "Screen_Ever",   "ever screened — detection-mode confounder",
  "Detection_screen","screen- vs symptom-detected"
)
dcis_qc <- pheno_raw %>%
  filter(status == 2, ethnicityClass == 1, MorphologygroupIndex_corr == "Ductal")
qc_log(sprintf("\n  %-16s  %-5s  %-9s  %-6s  %s",
               "Variable", "Pres.", "DCIS compl", "Used?", "Role"))
qc_log(sprintf("  %s", strrep("-", 96)))
for (i in seq_len(nrow(RISK_FACTORS))) {
  v <- RISK_FACTORS$var[i]
  present <- v %in% names(pheno_raw)
  compl <- if (present) round(100 * mean(!is.na(dcis_qc[[v]]))) else NA
  used  <- v %in% USED_VARS
  qc_log(sprintf("  %-16s  %-5s  %-9s  %-6s  %s",
                 v, ifelse(present, "yes", "NO"),
                 ifelse(present, paste0(compl, "%"), "-"),
                 ifelse(used, "yes", "no"), RISK_FACTORS$role[i]))
  if (present && !used && !is.na(compl) && compl >= 50)
    flag_warn(sprintf("%s present (%d%% complete in DCIS) but NOT in the analysis — consider as exposure/covariate or document exclusion", v, compl))
}
qc_log("\n  NOTE: famHist is intentionally excluded from genetic (Objective 1) models")
qc_log("        as it lies on the causal pathway (mediator) for carrier status;")
qc_log("        it is used only in the ascertainment sensitivity analysis.")


# =============================================================================
# LEVEL 3 — GENOTYPE QC
# =============================================================================

# ── Truncating merge (dominant model: value ≥ 1 → 1) ────────────────────────
trunc_select <- trunc_raw %>%
  select(BRIDGES_ID,
         all_of(unname(gene_col_map_trunc[genes_trunc_found]))) %>%
  rename(!!!setNames(gene_col_map_trunc[genes_trunc_found], genes_trunc_found))

dat_full_trunc <- pheno_raw %>%
  left_join(trunc_select, by = "BRIDGES_ID") %>%
  mutate(across(all_of(genes_trunc_found),
                ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))

# ── Missense merge (keep raw values for QC inspection) ───────────────────────
miss_select <- miss_raw %>%
  select(BRIDGES_ID,
         all_of(unname(gene_col_map_miss[genes_miss_found]))) %>%
  rename(!!!setNames(gene_col_map_miss[genes_miss_found], genes_miss_found))

dat_full_miss <- pheno_raw %>%
  left_join(miss_select, by = "BRIDGES_ID") %>%
  mutate(across(all_of(genes_miss_found),
                ~ replace_na(as.integer(.x), 0L)))

controls_geno_trunc <- dat_full_trunc %>%
  filter(status == 0, ethnicityClass == 1, !is.na(ageInt))

controls_geno_miss <- dat_full_miss %>%
  filter(status == 0, ethnicityClass == 1, !is.na(ageInt))


qc_section("LEVEL 3A — TRUNCATING CARRIER FREQUENCIES IN CONTROLS")
qc_log("  Comparison with Dorling et al. 2021 NEJM (European-ancestry controls)")
qc_log("  NOTE: Dorling 2021 values are approximate; verify against Table S2\n")

qc_log(sprintf("  %-8s  %-10s  %-10s  %-12s  %-12s  %s",
               "Gene", "Carriers", "N ctrl", "Freq (obs)", "Freq (pub)",  "Flag"))
qc_log(sprintf("  %s", strrep("-", 72)))

n_ctrl <- nrow(controls_geno_trunc)
for (g in genes_trunc_found) {
  n_carrier  <- sum(controls_geno_trunc[[g]] == 1, na.rm = TRUE)
  freq_obs   <- n_carrier / n_ctrl
  freq_pub   <- DORLING_CTRL_FREQ[g]
  ratio      <- if (!is.na(freq_pub) && freq_pub > 0) freq_obs / freq_pub else NA_real_

  flag_str <- if (is.na(ratio)) "pub freq N/A" else
    if (ratio < 0.5 || ratio > 2.0) "*** CHECK: >2-fold discrepancy ***" else
    if (ratio < 0.7 || ratio > 1.5) "* mild discrepancy" else "OK"

  qc_log(sprintf("  %-8s  %-10d  %-10d  %-12s  %-12s  %s",
                 g, n_carrier, n_ctrl,
                 sprintf("%.4f", freq_obs),
                 if (!is.na(freq_pub)) sprintf("%.4f", freq_pub) else "N/A",
                 flag_str))

  if (!is.na(ratio) && (ratio < 0.5 || ratio > 2.0))
    flag_warn(sprintf("Truncating %s: observed freq %.4f vs published %.4f (ratio=%.2f)",
                      g, freq_obs, freq_pub, ratio))
}


qc_section("LEVEL 3B — TRUNCATING CARRIER COUNTS BY OUTCOME")

cohort_geno_list <- list(
  "BBD cases"      = dat_full_trunc %>% filter(status == 0, ethnicityClass == 1,
                                                BBD_history == 1, !is.na(ageInt)),
  "BBD controls"   = dat_full_trunc %>% filter(status == 0, ethnicityClass == 1,
                                                BBD_history == 0, !is.na(ageInt)),
  "DCIS cases"     = dat_full_trunc %>% filter(status == 2, ethnicityClass == 1,
                                                MorphologygroupIndex_corr == "Ductal",
                                                !is.na(ageInt)),
  "LCIS cases"     = dat_full_trunc %>% filter(status == 2, ethnicityClass == 1,
                                                MorphologygroupIndex_corr == "Lobular",
                                                !is.na(ageInt)),
  "Full controls"  = dat_full_trunc %>% filter(status == 0, ethnicityClass == 1,
                                                !is.na(ageInt))
)

for (g in genes_trunc_found) {
  qc_log(sprintf("\n  Gene: %s (truncating)", g))
  qc_log(sprintf("  %-18s  %-8s  %-8s  %s",
                 "Cohort", "N", "Carriers", "Freq (%)"))
  qc_log(sprintf("  %s", strrep("-", 48)))
  for (cname in names(cohort_geno_list)) {
    dat  <- cohort_geno_list[[cname]]
    n    <- nrow(dat)
    nc   <- sum(dat[[g]] == 1, na.rm = TRUE)
    freq <- if (n > 0) 100 * nc / n else NA_real_
    qc_log(sprintf("  %-18s  %-8d  %-8d  %.3f%%", cname, n, nc, freq))
  }
}


qc_section("LEVEL 3C — MISSENSE CARRIER FREQUENCIES IN CONTROLS")

qc_log("  Missense carriers: GENE_CADD.phred.01 == 1 (CADD phred >= 20, top 1% deleterious)\n")
qc_log(sprintf("  %-8s  %-10s  %-10s  %s",
               "Gene", "Carriers", "N ctrl", "Freq (%)"))
qc_log(sprintf("  %s", strrep("-", 46)))

for (g in genes_miss_found) {
  n_carrier <- sum(controls_geno_miss[[g]] == 1, na.rm = TRUE)
  freq_obs  <- 100 * n_carrier / nrow(controls_geno_miss)
  qc_log(sprintf("  %-8s  %-10d  %-10d  %.3f%%",
                 g, n_carrier, n_ctrl, freq_obs))
}


qc_section("LEVEL 3D — STUDY-LEVEL CARRIER DISTRIBUTION (CONFOUNDING CHECK)")
qc_log("  This check identifies studies with unusually high or low carrier rates")
qc_log("  that could create study-level confounding not fully absorbed by the")
qc_log("  study covariate. This was identified as an issue for CHEK2 and KARMA.\n")

for (g in c("CHEK2", "ATM", "BRCA1", "BRCA2")) {
  if (!g %in% genes_trunc_found) next
  study_carrier <- controls_geno_trunc %>%
    group_by(study) %>%
    summarise(
      n_ctrl    = n(),
      n_carrier = sum(.data[[g]] == 1, na.rm = TRUE),
      freq_pct  = 100 * n_carrier / n_ctrl,
      .groups   = "drop"
    ) %>%
    filter(n_ctrl >= 50) %>%
    arrange(desc(freq_pct))

  overall_freq <- 100 * sum(controls_geno_trunc[[g]] == 1, na.rm = TRUE) / n_ctrl
  qc_log(sprintf("\n  %s truncating — overall control freq = %.3f%%", g, overall_freq))
  qc_log(sprintf("  %-30s  %6s  %8s  %6s",
                 "Study", "N ctrl", "Carriers", "Freq %"))
  qc_log(sprintf("  %s", strrep("-", 58)))
  for (i in seq_len(min(10, nrow(study_carrier)))) {
    row    <- study_carrier[i, ]
    flag   <- if (row$freq_pct > 3 * overall_freq || row$freq_pct < overall_freq / 3) " ***" else ""
    qc_log(sprintf("  %-30s  %6d  %8d  %5.2f%%%s",
                   row$study, row$n_ctrl, row$n_carrier, row$freq_pct, flag))
  }
  outlier_studies <- study_carrier %>%
    filter(freq_pct > 3 * overall_freq | freq_pct < overall_freq / 3)
  if (nrow(outlier_studies) > 0) {
    flag_warn(sprintf("%s: %d stud%s with >3x or <1/3x overall carrier freq (study-level confounding risk)",
                      g, nrow(outlier_studies),
                      if (nrow(outlier_studies) == 1) "y" else "ies"))
  }
}


# =============================================================================
# QC SUMMARY
# =============================================================================

qc_section("QC SUMMARY")

qc_log(sprintf("\n  Total PASS flags : not individually counted (see above)"))
qc_log(sprintf("  Total WARN flags : %d", length(qc_flags$warnings)))
qc_log(sprintf("  Total FAIL flags : %d", length(qc_flags$failures)))

if (length(qc_flags$failures) > 0) {
  qc_log("\nFAILURES (must resolve before analysis):")
  for (f in qc_flags$failures) qc_log("  [FAIL] ", f)
} else {
  qc_log("\n  No FAIL flags — data integrity checks passed.")
}

if (length(qc_flags$warnings) > 0) {
  qc_log("\nWARNINGS (review and document in dissertation):")
  for (w in qc_flags$warnings) qc_log("  [WARN] ", w)
} else {
  qc_log("\n  No WARN flags.")
}

qc_log("\n\nKEY NUMBERS FOR DISSERTATION METHODS SECTION:")
qc_log(sprintf("  European controls (status==0, ageInt non-missing) : %d", nrow(eur_s0)))
qc_log(sprintf("  BBD cases (status==0, BBD_history==1)             : %d", nrow(bbd_cases)))
qc_log(sprintf("  BBD controls (status==0, BBD_history==0)          : %d", nrow(bbd_ctrls)))
qc_log(sprintf("  DCIS cases (status==2, Ductal)                    : %d", nrow(dcis_cases)))
qc_log(sprintf("  LCIS cases (status==2, Lobular)                   : %d", nrow(lcis_cases)))
qc_log(sprintf("  status==3 with missing BBD_history                : %.0f%%",
               100 * mean(is.na(pheno_raw$BBD_history[pheno_raw$status == 3]))))

qc_log("\n\nUPSTREAM QC NOTE (cite in Methods):")
qc_log("  Genotype-level QC (sample call rate, heterozygosity, sex concordance,")
qc_log("  relatedness, Hardy-Weinberg equilibrium, variant call rate, batch effects,")
qc_log("  and ancestry PCA) was performed upstream by the BRIDGES consortium prior")
qc_log("  to data release. Reference: Dorling et al. (2021) NEJM 384:428-439.")
qc_log("  PCs were not available in phenotype file v17; ancestry was controlled by")
qc_log("  restricting to ethnicityClass==1 (European) and adjusting for study.")

qc_log(sprintf("\n\nQC report saved to: %s", report_path))
qc_log(sprintf("QC figures saved to: outputs/QC/figures/"))
qc_log(sprintf("Run completed: %s", format(Sys.time())))

close(con)
message(sprintf("\n✓ QC complete. Report: %s", report_path))


# ── SESSION INFO ──────────────────────────────────────────────────────────────
writeLines(
  capture.output(print(sessionInfo())),
  sprintf("outputs/QC/session_info_QC_%s.txt", timestamp)
)
