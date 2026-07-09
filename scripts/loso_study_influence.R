# =============================================================================
# Leave-One-Study-Out (LOSO) Sensitivity — DCIS truncating & missense
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil Population Health Sciences — 2025/26
#
# PURPOSE (addresses supervisor point #3: "study effect / is the result
#          dominated by a single contributing study?")
#   BRIDGES pools ~28 case-contributing studies. Adjusting for `study`
#   removes additive study effects, but a single large or atypical study
#   could still drive an association. This script re-estimates each DCIS
#   gene association 28 times, each time removing ONE entire contributing
#   study (its cases AND its controls), and reports:
#     - the full-data (primary) OR
#     - the min / max OR across all leave-one-out refits
#     - whether p < 0.05 holds in EVERY leave-one-out refit
#     - which single study's removal moves the OR the most
#
#   A finding is "robust to study influence" if the OR range is narrow and
#   p stays < 0.05 no matter which study is dropped. A finding is "study-
#   driven" if removing one study collapses (or creates) the association.
#
# Built to mirror sensitivity_familial_population.R exactly (same cohort
# definitions, same Firth estimator, same carrier construction).
#
# DATA GOVERNANCE: BCAC/BRIDGES approved secondary use; outputs are
#   non-disclosive summary statistics only.
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c("tidyverse", "logistf", "writexl")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({
  library(tidyverse)
  library(logistf)
  library(writexl)
})

# ── 1. PATHS AND CONSTANTS ────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
MISS_FILE  <- "concept_807_zhang_bridges_missense.csv"

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

GENES <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
           "BARD1", "RAD51C", "RAD51D", "TP53")

# The FDR-significant DCIS findings we most need to defend (for the summary).
PRIMARY_TARGETS <- tribble(
  ~vtype,        ~gene,
  "truncating",  "ATM",
  "truncating",  "CHEK2",
  "truncating",  "BRCA2",
  "missense",    "CHEK2",
  "missense",    "TP53"
)

# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
message("Loading phenotype file...")
pheno <- suppressWarnings(read_delim(
  PHENO_FILE, delim = "\t", show_col_types = FALSE,
  na = c("", "NA", "888", "777", "999")))

message("Loading truncating variants file...")
trunc <- read_csv(TRUNC_FILE, show_col_types = FALSE, na = c("", "NA")) %>%
  select(-any_of("...1"))

message("Loading missense variants file...")
miss <- read_csv(MISS_FILE, show_col_types = FALSE, na = c("", "NA")) %>%
  select(-any_of("...1"))

# Build a BRIDGES_ID + 9-gene 0/1 carrier table (identical to sensitivity script).
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

carriers_t <- build_carriers(trunc, "_truncating$",        "_t")
carriers_m <- build_carriers(miss,  "_CADD\\.phred\\.01$", "_m")

dat <- pheno %>%
  filter(ethnicityClass == 1, !is.na(ageInt)) %>%
  left_join(carriers_t, by = "BRIDGES_ID") %>%
  left_join(carriers_m, by = "BRIDGES_ID") %>%
  mutate(across(c(ends_with("_t"), ends_with("_m")), ~ replace_na(.x, 0L)))

# ── 3. COHORTS ────────────────────────────────────────────────────────────────
controls <- dat %>% filter(status == 0)
dcis     <- dat %>% filter(status == 2, MorphologygroupIndex_corr == "Ductal")

case_studies <- sort(unique(dcis$study))
message(sprintf("\nDCIS cases: %d across %d studies", nrow(dcis), length(case_studies)))

# ── 4. FIRTH FIT ON A GIVEN DATA SUBSET ───────────────────────────────────────
# drop_study = NULL -> full data (primary). Otherwise remove that study entirely
# (its cases AND controls). Returns OR/CI/p plus carrier counts after the drop.
fit_loso <- function(gene, vtype, drop_study = NULL) {
  carrier_col <- paste0(gene, if (vtype == "missense") "_m" else "_t")

  cs <- dcis; ct <- controls
  if (!is.null(drop_study)) {
    cs <- cs %>% filter(study != drop_study)
    ct <- ct %>% filter(study != drop_study)
  }

  d <- bind_rows(ct %>% mutate(outcome = 0L),
                 cs %>% mutate(outcome = 1L)) %>%
    mutate(carrier = .data[[carrier_col]], study = droplevels(factor(study)))

  n_cc <- sum(d$carrier == 1 & d$outcome == 1)
  n_ct <- sum(d$carrier == 1 & d$outcome == 0)

  base <- tibble(gene, vtype,
                 dropped = if (is.null(drop_study)) "(none — full data)" else drop_study,
                 n_case = sum(d$outcome == 1),
                 carr_case = n_cc, carr_ctrl = n_ct)

  if (n_cc == 0) {
    return(base %>% mutate(OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                           p_value = NA_real_, note = "0 carriers in cases"))
  }

  fit <- tryCatch(
    logistf(outcome ~ carrier + ageInt + study, data = d, firth = TRUE, pl = TRUE, plconf = 2),
    error = function(e) NULL)
  if (is.null(fit)) {
    return(base %>% mutate(OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                           p_value = NA_real_, note = "fit failed"))
  }
  i <- which(names(coef(fit)) == "carrier")
  base %>% mutate(OR = exp(coef(fit)[i]),
                  CI_low = exp(fit$ci.lower[i]), CI_high = exp(fit$ci.upper[i]),
                  p_value = fit$prob[i], note = "")
}

# ── 5. RUN: full fit + drop each case-contributing study, per gene × vtype ─────
# CHECKPOINT: the ~520 Firth refits take ~12 min. Save immediately after fitting
# (and reload on re-run) so the summary/export logic can be iterated for free.
vtypes <- c("truncating", "missense")
CHECKPOINT <- "outputs/tables/loso_all_fits.rds"

if (file.exists(CHECKPOINT)) {
  message("\nLoading cached LOSO fits from ", CHECKPOINT,
          " (delete this file to force a refit).")
  all_loso <- readRDS(CHECKPOINT)
} else {
  message("\nRunning LOSO (full + ", length(case_studies),
          " single-study drops) for 9 genes × 2 variant classes...")
  all_loso <- map_dfr(vtypes, function(vt) {
    map_dfr(GENES, function(g) {
      full <- fit_loso(g, vt, NULL)
      drops <- map_dfr(case_studies, function(s) fit_loso(g, vt, s))
      bind_rows(full, drops)
    })
  })
  saveRDS(all_loso, CHECKPOINT)
  message("Checkpoint saved: ", CHECKPOINT)
}

# ── 6. PER-GENE ROBUSTNESS SUMMARY ────────────────────────────────────────────
# Robust to genes where every leave-one-out fit is NA (e.g. RAD51D truncating,
# 0 carriers among DCIS cases): such genes return NA across the summary columns
# rather than crashing on an empty which.max().
summarise_group <- function(df) {
  full <- df %>% filter(dropped == "(none — full data)")
  loso <- df %>% filter(dropped != "(none — full data)")
  valid <- loso %>% filter(!is.na(OR))
  OR_full <- full$OR[1]; p_full <- full$p_value[1]

  if (nrow(valid) == 0) {
    return(tibble(OR_full, p_full,
                  OR_min = NA_real_, OR_max = NA_real_, p_max = NA_real_,
                  n_loso_fits = 0L, p_lt05_always = NA,
                  most_infl_study = NA_character_, most_infl_OR = NA_real_))
  }
  k <- which.max(abs(valid$OR - OR_full))
  tibble(
    OR_full, p_full,
    OR_min = min(valid$OR), OR_max = max(valid$OR),
    p_max  = max(valid$p_value, na.rm = TRUE),
    n_loso_fits = nrow(valid),
    p_lt05_always = all(valid$p_value < 0.05),
    most_infl_study = valid$dropped[k],
    most_infl_OR    = valid$OR[k]
  )
}

summ <- all_loso %>%
  group_by(gene, vtype) %>%
  group_modify(~ summarise_group(.x)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, levels = GENES),
         vtype = factor(vtype, levels = vtypes)) %>%
  arrange(vtype, gene)

cat("\n", strrep("=", 90), "\n",
    "DCIS LEAVE-ONE-STUDY-OUT SUMMARY (each row = gene; range over all single-study drops)\n",
    strrep("=", 90), "\n", sep = "")
summ %>%
  mutate(`OR full (p)`    = sprintf("%.2f (p=%.3f)", OR_full, p_full),
         `OR range (LOSO)`= sprintf("%.2f – %.2f", OR_min, OR_max),
         `p<0.05 always?` = if_else(p_lt05_always, "YES", "NO"),
         `worst p`        = sprintf("%.3f", p_max),
         `most influential drop` = sprintf("%s -> OR %.2f", most_infl_study, most_infl_OR)) %>%
  select(vtype, gene, `OR full (p)`, `OR range (LOSO)`,
         `p<0.05 always?`, `worst p`, `most influential drop`) %>%
  as.data.frame() %>% print(right = FALSE)

cat("\n--- FOCUS: the 5 FDR-significant DCIS findings ---\n")
PRIMARY_TARGETS %>%
  left_join(summ, by = c("vtype", "gene")) %>%
  mutate(verdict = if_else(p_lt05_always & OR_min > 1,
                           "ROBUST (p<0.05 in all drops, OR>1 throughout)",
                           "CHECK (significance or direction lost in >=1 drop)")) %>%
  select(vtype, gene, OR_full, OR_min, OR_max, p_max, most_infl_study, verdict) %>%
  as.data.frame() %>% print(right = FALSE)

# ── 7. EXPORT ─────────────────────────────────────────────────────────────────
out_path <- "outputs/tables/loso_study_influence.xlsx"
write_xlsx(
  list(
    "LOSO_summary" = summ %>% as.data.frame(),
    "LOSO_full_long" = all_loso %>% as.data.frame()
  ),
  path = out_path
)
message(sprintf("\nExported: %s", out_path))

sink("outputs/session_info/session_info_loso.txt")
cat("LOSO Study-Influence Sensitivity — Session Info\n")
cat(format(Sys.time()), "\n\n"); sessionInfo()
sink()
message("Done.")
