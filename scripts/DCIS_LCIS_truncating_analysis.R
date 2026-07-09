# =============================================================================
# Dissertation Analysis — DCIS and LCIS Truncating Variant Analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil — Objective 1 (DCIS/LCIS component)
#
# DATA GOVERNANCE
#   Data accessed under BCAC/BRIDGES approved secondary use agreement.
#   All analysis conducted on non-identifiable summary-level outputs only.
#
# GENOTYPE QC PROVENANCE
#   Variant-level QC applied upstream by BRIDGES consortium prior to release.
#   Reference: Dorling et al. (2021) NEJM 384:428-439.
#   No additional variant-level QC applied here.
#   PCs not available in v17; ancestry controlled by ethnicityClass==1 + study.
#
# COHORT DEFINITIONS:
#   DCIS cases    : status==2, MorphologygroupIndex_corr=="Ductal",
#                   ethnicityClass==1, !is.na(ageInt)
#   LCIS cases    : status==2, MorphologygroupIndex_corr=="Lobular",
#                   ethnicityClass==1, !is.na(ageInt)
#   Controls      : status==0, ethnicityClass==1, !is.na(ageInt)
#   (same control pool as missense_analysis_corrected.R)
#
# MODEL: Firth's penalised logistic regression (logistf, pl=FALSE → Wald CIs;
#        profile-likelihood hangs on the 0-carrier LCIS genes — see section 7.
#        Confirmed materially identical to PL for the significant genes.)
# COVARIATES: ageInt + study (categorical)
# FDR: BH across 9 genes, separately for DCIS and LCIS
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c("tidyverse", "logistf", "writexl", "ggplot2", "scales")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(tidyverse)
library(logistf)
library(writexl)
library(ggplot2)
library(scales)

# ── 1. FILE PATHS ─────────────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"

# Ensure output directories exist before writing
dir.create("outputs/tables",       recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures",      recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/session_info", recursive = TRUE, showWarnings = FALSE)

# Working-directory guard: scripts use paths relative to the project root.
for (f in c(PHENO_FILE, TRUNC_FILE)) {
  if (!file.exists(f))
    stop("Cannot find '", f, "'. Run this script from the project root ",
         "(the directory containing the BRIDGES data files).", call. = FALSE)
}

# ── 2. GENE LIST ──────────────────────────────────────────────────────────────
GENES <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
           "BARD1", "RAD51C", "RAD51D", "TP53")

# ── 3. LOAD DATA ──────────────────────────────────────────────────────────────
message("Loading phenotype file...")
pheno_raw <- read_delim(PHENO_FILE, delim = "\t", show_col_types = FALSE,
                        na = c("", "NA", "888", "777", "999"))

message("Loading truncating variants file...")
# Genotype files are expected to contain only 0 (non-carrier) and 1 (carrier).
# NA codes 888/777/999 are NOT added here as they are not expected in genotype
# columns; blank/NA strings are handled by na = c("", "NA").
trunc_raw <- read_csv(TRUNC_FILE, show_col_types = FALSE,
                      na = c("", "NA")) %>%
  select(-any_of("...1"))

# ── 4. AUTO-DETECT GENE COLUMNS ───────────────────────────────────────────────
detect_gene_col <- function(gene, col_names) {
  exact  <- col_names[tolower(col_names) == tolower(gene)]
  if (length(exact) > 0) return(exact[1])
  prefix <- col_names[grepl(paste0("^", gene, "(_|$)"), col_names, ignore.case = TRUE)]
  if (length(prefix) > 0) return(prefix[1])
  NA_character_
}

gene_col_map <- setNames(
  sapply(GENES, detect_gene_col, col_names = names(trunc_raw)), GENES)
genes_found <- names(gene_col_map)[!is.na(gene_col_map)]
message("Gene columns mapped: ", paste(genes_found, collapse = ", "))

# ── 5. MERGE ──────────────────────────────────────────────────────────────────
trunc_select <- trunc_raw %>%
  select(BRIDGES_ID, all_of(unname(gene_col_map[genes_found]))) %>%
  rename(!!!setNames(gene_col_map[genes_found], genes_found))

# Join-coverage guard: every phenotype participant must have a genotype record,
# else left_join + replace_na(0) would silently code them as non-carriers.
n_unmatched <- pheno_raw %>%
  anti_join(trunc_select, by = "BRIDGES_ID") %>%
  nrow()
if (n_unmatched > 0)
  stop(sprintf("%d phenotype IDs have no truncating-genotype record and would be ",
               n_unmatched),
       "silently coded as non-carriers. Aborting.", call. = FALSE)

dat_full <- pheno_raw %>%
  left_join(trunc_select, by = "BRIDGES_ID") %>%
  # Dominant model: any carrier (≥ 1 truncating variant) → 1, non-carrier → 0
  mutate(across(all_of(genes_found),
                ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))

message(sprintf("Merged dataset: %d rows", nrow(dat_full)))

# ── 6. BUILD COHORTS ──────────────────────────────────────────────────────────
controls <- dat_full %>%
  filter(status == 0, ethnicityClass == 1, !is.na(ageInt)) %>%
  mutate(outcome = 0L)

dcis_cases <- dat_full %>%
  filter(status == 2, ethnicityClass == 1,
         MorphologygroupIndex_corr == "Ductal", !is.na(ageInt)) %>%
  mutate(outcome = 1L)

lcis_cases <- dat_full %>%
  filter(status == 2, ethnicityClass == 1,
         MorphologygroupIndex_corr == "Lobular", !is.na(ageInt)) %>%
  mutate(outcome = 1L)

dat_dcis <- bind_rows(controls, dcis_cases)
dat_lcis <- bind_rows(controls, lcis_cases)

message(sprintf("Controls     : %d  [expected 34,424]", nrow(controls)))
message(sprintf("DCIS cases   : %d  [expected 1,663]",  nrow(dcis_cases)))
message(sprintf("LCIS cases   : %d  [expected 139]",    nrow(lcis_cases)))

# ── 7. FIRTH'S REGRESSION HELPER ──────────────────────────────────────────────
run_firth <- function(data, gene) {
  n_carrier_case <- sum(data[[gene]] == 1 & data$outcome == 1, na.rm = TRUE)
  n_carrier_ctrl <- sum(data[[gene]] == 1 & data$outcome == 0, na.rm = TRUE)
  n_case         <- sum(data$outcome == 1, na.rm = TRUE)
  n_ctrl         <- sum(data$outcome == 0, na.rm = TRUE)

  # ── Pre-check: skip if no carriers in either group ────────────────────────
  # Prevents a LAPACK error when the predictor is entirely absent (e.g. TP53
  # with 0 carriers in both LCIS cases and controls). Reported as NA.
  if ((n_carrier_case + n_carrier_ctrl) == 0) {
    message(sprintf("  %s: 0 carriers in both groups — skipped (reported as NA)", gene))
    return(tibble(gene, n_case, n_ctrl,
                  n_carrier_case = 0L, n_carrier_ctrl = 0L,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "no carriers"))
  }

  model_data <- data %>%
    select(outcome, carrier = all_of(gene), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()

  # ── CI method: Wald (pl = FALSE) ──────────────────────────────────────────
  # Wald CIs are used here because several LCIS genes have 0 carriers among
  # cases (e.g. BARD1, RAD51C, TP53); profile-likelihood (pl = TRUE) profiling
  # of an almost-unidentified coefficient on these genes is extremely slow and
  # can hang. For the genes that drive the findings (ATM, CHEK2, BRCA2) Wald and
  # profile-likelihood CIs are materially identical because carriers are ample;
  # this was confirmed in a targeted pl = TRUE check (see compare_ci_methods.R).
  # The very sparse genes are not interpretable under either CI method.
  fit <- tryCatch(
    logistf(outcome ~ carrier + ageInt + study,
            data = model_data, firth = TRUE, pl = FALSE),
    error = function(e) {
      message(sprintf("  Firth failed for %s: %s", gene, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) {
    return(tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "failed"))
  }

  carrier_idx <- which(names(coef(fit)) == "carrier")
  tibble(gene,
         n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
         OR      = exp(coef(fit)[carrier_idx]),
         CI_low  = exp(fit$ci.lower[carrier_idx]),
         CI_high = exp(fit$ci.upper[carrier_idx]),
         p_value = fit$prob[carrier_idx],
         method  = "Firth")
}

# ── 8. RUN ANALYSES ───────────────────────────────────────────────────────────
run_panel <- function(data, label) {
  message(sprintf("\nRunning: %s", label))
  map_dfr(genes_found, ~ run_firth(data, .x)) %>%
    mutate(analysis = label,
           p_fdr    = p.adjust(p_value, method = "BH"),
           fdr_sig  = p_fdr < 0.05) %>%
    select(analysis, gene, n_case, n_ctrl,
           n_carrier_case, n_carrier_ctrl,
           OR, CI_low, CI_high, p_value, p_fdr, fdr_sig, method)
}

res_dcis <- run_panel(dat_dcis, "DCIS vs Controls (truncating)")
res_lcis <- run_panel(dat_lcis, "LCIS vs Controls (truncating)")

# ── 9. FORMAT TABLE ───────────────────────────────────────────────────────────
fmt_table <- function(df) {
  df %>%
    mutate(
      `OR (95% CI)` = case_when(
        is.na(OR) ~ "—",
        TRUE      ~ sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)
      ),
      `p-value` = case_when(
        is.na(p_value)  ~ "—",
        p_value < 0.001 ~ "<0.001",
        TRUE            ~ sprintf("%.3f", p_value)
      ),
      `p-FDR (BH)` = case_when(
        is.na(p_fdr)  ~ "—",
        p_fdr < 0.001 ~ "<0.001",
        TRUE          ~ sprintf("%.3f", p_fdr)
      ),
      `FDR sig.` = if_else(fdr_sig, "Yes*", "No", missing = "—")
    ) %>%
    select(Gene = gene,
           `N (cases)`           = n_case,
           `N (controls)`        = n_ctrl,
           `Carriers (cases)`    = n_carrier_case,
           `Carriers (controls)` = n_carrier_ctrl,
           `OR (95% CI)`,
           `p-value`,
           `p-FDR (BH)`,
           `FDR sig.`)
}

# Print to console
for (lbl in c("DCIS vs Controls (truncating)", "LCIS vs Controls (truncating)")) {
  cat(strrep("─", 70), "\n", lbl, "\n", strrep("─", 70), "\n")
  df <- if (grepl("DCIS", lbl)) res_dcis else res_lcis
  print(fmt_table(df), n = Inf)
  cat("\n")
}

# ── 10. EXPORT EXCEL ──────────────────────────────────────────────────────────
write_xlsx(
  list(
    "DCIS_truncating" = fmt_table(res_dcis),
    "LCIS_truncating" = fmt_table(res_lcis),
    "Raw_results"     = bind_rows(res_dcis, res_lcis)
  ),
  path = "outputs/tables/DCIS_LCIS_truncating_results.xlsx"
)
message("Excel exported: outputs/tables/DCIS_LCIS_truncating_results.xlsx")

# ── 11. FOREST PLOT FUNCTION ──────────────────────────────────────────────────
# NOTE: the canonical figures for the dissertation are now regenerated centrally
# — in one unified house style — by scripts/build_exhibits.R, which calls the
# shared forest_plot() helper (scripts/_exhibits_helpers.R) on the saved result
# tables. This inline version is retained only for standalone runs of this script.
make_forest <- function(res_df, title_text, filename) {
  plot_df <- res_df %>%
    mutate(
      gene_f       = factor(gene, levels = rev(GENES)),
      colour_group = case_when(
        is.na(OR) ~ "Failed",
        fdr_sig   ~ "FDR significant",
        TRUE      ~ "Not significant"
      ),
      or_label      = case_when(
        is.na(OR) ~ "Model failed",
        TRUE      ~ sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)
      ),
      carrier_label = sprintf("%d / %d", n_carrier_case, n_carrier_ctrl)
    )

  x_min <- min(c(plot_df$CI_low,  0.3), na.rm = TRUE) * 0.75
  # Cap the panel max so artifact CIs (e.g. LCIS TP53/BARD1, upper bound > 200)
  # do not compress the log axis and shove the side text off the right edge.
  x_max <- min(max(c(plot_df$CI_high, 5.0), na.rm = TRUE), 25) * 1.20
  # Clip the drawn error bar to the data region (full CI is still printed in the
  # OR text column); this keeps over-cap artifact bars from running into the text.
  plot_df <- plot_df %>%
    mutate(CI_high_draw = pmin(CI_high, x_max),
           clipped      = !is.na(CI_high) & CI_high > x_max)

  colour_vals <- c("FDR significant" = "#B03A2E",
                   "Not significant" = "#1A5276",
                   "Failed"          = "grey70")
  shape_vals  <- c("FDR significant" = 18,
                   "Not significant" = 18,
                   "Failed"          = 1)

  x_carrier <- x_max * 1.30
  x_or      <- x_max * 2.40

  p <- ggplot(plot_df, aes(y = gene_f)) +
    geom_vline(xintercept = 1, linetype = "longdash",
               colour = "grey55", linewidth = 0.5) +
    geom_errorbar(
      aes(xmin = CI_low, xmax = CI_high_draw, colour = colour_group),
      width = 0.28, linewidth = 0.8, na.rm = TRUE, orientation = "y"
    ) +
    # small arrow where the CI was clipped at the panel edge (artifact genes)
    geom_text(data = ~ filter(.x, clipped),
              aes(x = x_max, label = "→"),
              hjust = 0, vjust = 0.35, size = 3.2, colour = "#1A5276", na.rm = TRUE) +
    geom_point(
      aes(x = OR, colour = colour_group, shape = colour_group),
      size = 4, na.rm = TRUE
    ) +
    geom_text(aes(x = x_carrier, label = carrier_label),
              hjust = 0.5, size = 3, colour = "grey25") +
    geom_text(aes(x = x_or, label = or_label),
              hjust = 0, size = 3, colour = "grey15") +
    annotate("text", x = x_carrier, y = length(GENES) + 0.7,
             label = "Cases / Ctrl", size = 3, fontface = "bold", hjust = 0.5) +
    annotate("text", x = x_or, y = length(GENES) + 0.7,
             label = "OR (95% CI)", size = 3, fontface = "bold", hjust = 0) +
    scale_x_log10(
      breaks = c(0.3, 0.5, 1, 2, 5, 10, 20),
      labels = c("0.3", "0.5", "1", "2", "5", "10", "20"),
      limits = c(x_min, x_or * 1.9)
    ) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
    scale_colour_manual(values = colour_vals, name = NULL) +
    scale_shape_manual(values = shape_vals, name = NULL) +
    labs(title   = title_text,
         x       = "Odds Ratio (log scale, 95% Wald CI)",
         y       = NULL,
         caption = "Firth's penalised logistic regression (Wald CIs). Adjusted for age and study.\nBH-FDR across 9 genes. Red = FDR q < 0.05.") +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", size = 13),
      plot.caption     = element_text(size = 9, colour = "grey40")
    )

  ggsave(filename, plot = p, width = 12, height = 6.5, dpi = 300, bg = "white")
  message(sprintf("Forest plot saved: %s", filename))
  p
}

make_forest(res_dcis,
            "Truncating Variant Associations with DCIS",
            "outputs/figures/forest_DCIS_truncating.png")

make_forest(res_lcis,
            "Truncating Variant Associations with LCIS",
            "outputs/figures/forest_LCIS_truncating.png")

# ── 12. SESSION INFO ──────────────────────────────────────────────────────────
sink("outputs/session_info/session_info_DCIS_LCIS_truncating.txt")
cat("DCIS/LCIS Truncating Analysis — Session Info\n")
cat(format(Sys.time()), "\n\n")
print(sessionInfo())  # print() required inside sink() to write output to file
sink()
message("Session info saved.")
message("\nDone. All outputs in outputs/tables/ and outputs/figures/")
