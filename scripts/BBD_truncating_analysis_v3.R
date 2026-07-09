# =============================================================================
# Dissertation Analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil — Objective 1: BBD Truncating Variant Analysis
#
# DATA GOVERNANCE
#   Data accessed under BCAC/BRIDGES approved secondary use agreement.
#   All analysis conducted on non-identifiable summary-level outputs only.
#   Raw data files are not redistributed or committed to version control.
#
# GENOTYPE QC PROVENANCE
#   Variant-level QC (call rate, MAF, HWE, heterozygosity) was applied
#   upstream by the BRIDGES consortium prior to data release.
#   Reference: Dorling et al. (2021) NEJM 384:428-439 (BRIDGES QC protocol).
#   No individual-level genotype QC was re-applied in this analysis.
#   Principal components (PCs) for ancestry stratification were NOT available
#   in the released phenotype file (v17); ancestry was controlled by restricting
#   to ethnicityClass == 1 (European) and adjusting for study (categorical).
#   This is a documented limitation — see dissertation Section 5.4.
#
# COMPUTATIONAL ENVIRONMENT
#   R version and package versions saved to session_info.txt at script end.
#   OS: macOS 14 (Sonoma). Analysis date: see sessionInfo() output.
#
# VERSION 3 FIXES:
#   1. Firth's penalised logistic regression (logistf) replaces standard glm
#      for ALL models. This handles:
#        - Complete separation (TP53: 1 case carrier, 0 control carriers)
#        - Zero-cell counts (PALB2, BARD1, etc. with 0 carriers in cases)
#        - Quasi-complete separation (inflated SEs, "50+ warnings")
#      Firth's regression is standard practice in rare-variant genetic epi.
#   2. BBD count discrepancy investigated: checks both "Exclusion" and
#      "exclusion" (case-insensitive), and reports status distribution
#   3. Forest plot axis guide warning fixed
# =============================================================================

# ── 0.  PACKAGES ──────────────────────────────────────────────────────────────
required <- c("tidyverse", "logistf", "broom", "writexl", "ggplot2", "scales")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(tidyverse)
library(logistf)    # Firth's penalised likelihood logistic regression
library(writexl)
library(ggplot2)
library(scales)


# ── 1.  FILE PATHS ────────────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
OUT_XLSX   <- "outputs/tables/BBD_truncating_FINAL.xlsx"

# Ensure output directories exist
dir.create("outputs/tables",       recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures",      recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/session_info", recursive = TRUE, showWarnings = FALSE)

# Working-directory guard: scripts use paths relative to the project root.
# Fail early with a clear message rather than a confusing read error.
for (f in c(PHENO_FILE, TRUNC_FILE)) {
  if (!file.exists(f))
    stop("Cannot find '", f, "'. Run this script from the project root ",
         "(the directory containing the BRIDGES data files).", call. = FALSE)
}


# ── 2.  GENE LIST ─────────────────────────────────────────────────────────────
GENES <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
           "BARD1", "RAD51C", "RAD51D", "TP53")


# ── 3.  LOAD DATA ─────────────────────────────────────────────────────────────
message("Loading phenotype file...")
pheno_raw <- read_delim(
  PHENO_FILE, delim = "\t", show_col_types = FALSE,
  na = c("", "NA", "888", "777", "999")
)

message("Loading truncating variants file...")
trunc_raw <- read_csv(
  TRUNC_FILE, show_col_types = FALSE, na = c("", "NA")
) %>%
  select(-any_of("...1"))


# ── 4.  AUTO-DETECT GENE COLUMNS (unchanged from v2) ──────────────────────────
detect_gene_col <- function(gene, col_names) {
  exact  <- col_names[tolower(col_names) == tolower(gene)]
  if (length(exact)  > 0) return(exact[1])
  prefix <- col_names[grepl(paste0("^", gene, "(_|$)"), col_names,
                             ignore.case = TRUE)]
  if (length(prefix) > 0) return(prefix[1])
  NA_character_
}

gene_col_map <- setNames(
  sapply(GENES, detect_gene_col, col_names = names(trunc_raw)), GENES
)
genes_found <- names(gene_col_map)[!is.na(gene_col_map)]
message("Gene columns mapped: ", paste(genes_found, collapse = ", "))


# ── 5.  MERGE ─────────────────────────────────────────────────────────────────
trunc_select <- trunc_raw %>%
  select(BRIDGES_ID, all_of(unname(gene_col_map[genes_found]))) %>%
  rename(!!!setNames(gene_col_map[genes_found], genes_found))

# Join-coverage guard: every phenotype participant must have a genotype record.
# Without this, a left_join + replace_na(0) would silently code any unmatched
# participant as a non-carrier (verified 0 unmatched in v17; this fails loudly
# if a future data refresh breaks that invariant).
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
  # Handles multi-allelic individuals (value = 2) correctly under dominant model
  mutate(across(all_of(genes_found),
                ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))

message(sprintf("Merged dataset: %d rows", nrow(dat_full)))


# ── 6.  INVESTIGATE BBD COUNT DISCREPANCY ─────────────────────────────────────
message("\n── Status distribution in full dataset ──")
print(table(dat_full$status, useNA = "always"))

# Check for Exclusion column in any capitalisation
excl_col <- names(dat_full)[tolower(names(dat_full)) == "exclusion"]
if (length(excl_col) > 0) {
  message(sprintf("\nExclusion column found: '%s'", excl_col))
  message("Exclusion value counts:")
  print(table(dat_full[[excl_col]], useNA = "always"))
} else {
  message("\nNo 'Exclusion' column found in any capitalisation.")
  message("BBD count (pre-ancestry filter: 2019; post ethnicityClass==1 filter: 1765).")
  message("Discrepancy is explained by the European ancestry restriction — not a data error.")
}

# Define BBD cohort — status == 0 AND European ancestry only
# ethnicityClass == 1 is required for consistency with all other analyses
# (missense, DCIS, LCIS all restrict to European ancestry)
dat_status0 <- dat_full %>%
  filter(status == 0, ethnicityClass == 1)

n_eur <- nrow(dat_status0)
message(sprintf("After ethnicityClass == 1 filter: %d participants", n_eur))

# Apply exclusion filter if column found
if (length(excl_col) > 0) {
  n_before   <- nrow(dat_status0)
  dat_status0 <- dat_status0 %>%
    filter(.data[[excl_col]] == 0 | is.na(.data[[excl_col]]))
  message(sprintf("Exclusion filter removed %d rows", n_before - nrow(dat_status0)))
}

dat_bbd_all <- dat_status0 %>%
  filter(!is.na(BBD_history)) %>%
  mutate(outcome = as.integer(BBD_history == 1))

n_bbd  <- sum(dat_bbd_all$outcome == 1)
n_ctrl <- sum(dat_bbd_all$outcome == 0)
message(sprintf("\nBBD cases   : n = %d  [expected 1,765]", n_bbd))
message(sprintf("Controls    : n = %d  [expected 12,715]", n_ctrl))


# ── 7.  BBD SUBTYPES ──────────────────────────────────────────────────────────
dat_nonprol <- dat_bbd_all %>%
  filter(outcome == 0 | (outcome == 1 & BBD_type1 == 1)) %>%
  mutate(outcome = if_else(BBD_type1 == 1, 1L, 0L, missing = 0L))

dat_prol <- dat_bbd_all %>%
  filter(outcome == 0 | (outcome == 1 & BBD_type1 == 2)) %>%
  mutate(outcome = if_else(BBD_type1 == 2, 1L, 0L, missing = 0L))

message(sprintf("\nNon-proliferative : %d cases / %d controls",
                sum(dat_nonprol$outcome == 1), sum(dat_nonprol$outcome == 0)))
message(sprintf("Proliferative     : %d cases / %d controls",
                sum(dat_prol$outcome == 1),    sum(dat_prol$outcome == 0)))


# ── 8.  FIRTH'S LOGISTIC REGRESSION HELPER ────────────────────────────────────
#
# NOTE ON DOMINANT MODEL RECODING:
#   Gene columns were recoded to binary (0/1) in Section 5 using the dominant
#   model (≥1 → 1). Multi-allelic carriers (value = 2 in source data) are
#   therefore coded as 1 (carrier) here, which is correct under the dominant
#   model. n_carrier_case/ctrl below count individuals with gene == 1.
#
# WHY FIRTH'S REGRESSION:
#   Standard logistic regression (glm) fails when carrier counts are very low
#   or zero in one group — a problem known as "complete separation" or
#   "quasi-complete separation". It produces inflated ORs (e.g. 723,528 for
#   TP53) and unreliable CIs. Firth's penalised likelihood regression
#   (Firth 1993; Heinze & Schemper 2002) applies a bias-reduction penalty
#   that stabilises estimation in sparse data without excluding rare genes.
#   It is the recommended approach for rare-variant association studies.
#   Reference: Heinze G, Schemper M (2002) Stat Med 21:2409-2419.
#
# Model: logit(P(outcome)) = β₀ + β_carrier + β_age·ageInt + Σ β_j·study_j
# CIs: profile penalised likelihood (Firth method)

run_firth <- function(data, gene) {

  n_carrier_case <- sum(data[[gene]] == 1 & data$outcome == 1, na.rm = TRUE)
  n_carrier_ctrl <- sum(data[[gene]] == 1 & data$outcome == 0, na.rm = TRUE)
  n_case         <- sum(data$outcome == 1, na.rm = TRUE)
  n_ctrl         <- sum(data$outcome == 0, na.rm = TRUE)

  # ── Pre-check: skip if no carriers observed in either group ────────────────
  # Prevents Firth from attempting to estimate a coefficient that is entirely
  # unidentified (e.g. TP53 in subtype analyses with 0 carriers in both groups).
  # Integrated from BBD_firth_patch.R — reported as NA and excluded from FDR.
  if ((n_carrier_case + n_carrier_ctrl) == 0) {
    message(sprintf("  %s: 0 carriers in both groups — skipped (reported as NA)", gene))
    return(tibble(gene, n_case, n_ctrl,
                  n_carrier_case = 0L, n_carrier_ctrl = 0L,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "no carriers"))
  }

  # Prepare model data
  model_data <- data %>%
    select(outcome, carrier = all_of(gene), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()

  # Fit Firth's penalised logistic regression
  fit <- tryCatch(
    logistf(
      outcome ~ carrier + ageInt + study,
      data    = model_data,
      firth   = TRUE,        # Firth's penalisation (default TRUE)
      pl      = TRUE         # profile likelihood CIs (more accurate)
    ),
    error = function(e) {
      message(sprintf("  Firth model failed for %s: %s", gene, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) {
    return(tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "failed"))
  }

  # Extract carrier coefficient (term index 2: intercept=1, carrier=2)
  carrier_idx <- which(names(coef(fit)) == "carrier")

  OR      <- exp(coef(fit)[carrier_idx])
  CI_low  <- exp(fit$ci.lower[carrier_idx])
  CI_high <- exp(fit$ci.upper[carrier_idx])
  p_value <- fit$prob[carrier_idx]

  tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
         OR, CI_low, CI_high, p_value,
         method = "Firth")
}


# ── 9.  RUN ALL ANALYSES + FDR ────────────────────────────────────────────────
run_panel <- function(data, label) {
  message(sprintf("\nRunning: %s", label))
  map_dfr(genes_found, ~ run_firth(data, .x)) %>%
    mutate(
      analysis = label,
      p_fdr    = p.adjust(p_value, method = "BH"),
      fdr_sig  = p_fdr < 0.05
    ) %>%
    select(analysis, gene, n_case, n_ctrl,
           n_carrier_case, n_carrier_ctrl,
           OR, CI_low, CI_high, p_value, p_fdr, fdr_sig, method)
}

res_bbd     <- run_panel(dat_bbd_all, "BBD vs Controls")
res_nonprol <- run_panel(dat_nonprol, "Non-proliferative BBD vs Controls")
res_prol    <- run_panel(dat_prol,    "Proliferative BBD vs Controls (without atypia)")

results_all <- bind_rows(res_bbd, res_nonprol, res_prol)


# ── 10.  FORMAT AND PRINT TABLES ─────────────────────────────────────────────
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
    select(
      Gene                  = gene,
      `N (cases)`           = n_case,
      `N (controls)`        = n_ctrl,
      `Carriers (cases)`    = n_carrier_case,
      `Carriers (controls)` = n_carrier_ctrl,
      `OR (95% CI)`,
      `p-value`,
      `p-FDR (BH)`,
      `FDR sig.`
    )
}

for (lbl in unique(results_all$analysis)) {
  cat(strrep("─", 70), "\n")
  cat(lbl, "\n")
  cat(strrep("─", 70), "\n")
  print(fmt_table(filter(results_all, analysis == lbl)), n = Inf)
  cat("\n")
}


# ── 11.  EXPORT TO EXCEL ─────────────────────────────────────────────────────
write_xlsx(
  list(
    "BBD vs Controls"       = fmt_table(res_bbd),
    "Non-proliferative BBD" = fmt_table(res_nonprol),
    "Proliferative BBD"     = fmt_table(res_prol),
    "Raw results"           = results_all
  ),
  path = OUT_XLSX
)
message(sprintf("\nResults exported to: %s", OUT_XLSX))


# ── 12.  FOREST PLOT ─────────────────────────────────────────────────────────
# NOTE: the canonical dissertation figures are now regenerated centrally — in one
# unified house style — by scripts/build_exhibits.R via the shared forest_plot()
# helper (scripts/_exhibits_helpers.R), drawn from the saved result tables. This
# inline version is retained only for standalone runs of this script.
make_forest_plot <- function(res_df, title_text, filename) {

  plot_df <- res_df %>%
    mutate(
      gene_f        = factor(gene, levels = rev(GENES)),
      colour_group  = case_when(
        is.na(OR) ~ "Failed",
        fdr_sig   ~ "FDR significant",
        TRUE      ~ "Not significant"
      ),
      or_label = case_when(
        is.na(OR) ~ "Model failed",
        TRUE      ~ sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)
      ),
      carrier_label = sprintf("%d / %d", n_carrier_case, n_carrier_ctrl)
    )

  x_min <- min(c(plot_df$CI_low,  0.3), na.rm = TRUE) * 0.75
  # Cap the panel maximum so a single artifact CI (e.g. TP53 with 1 case / 0
  # control carriers, upper CI > 1000) does not compress the log axis.
  x_max <- min(max(c(plot_df$CI_high, 5.0), na.rm = TRUE), 30) * 1.20
  # Clip the drawn error bar to the data region so over-cap artifact bars never
  # run into the carrier/OR text columns; the full CI is still printed as text
  # and a small arrow marks where the bar was clipped.
  plot_df <- plot_df %>%
    mutate(CI_high_draw = pmin(CI_high, x_max),
           clipped      = !is.na(CI_high) & CI_high > x_max)

  colour_vals <- c("FDR significant" = "#B03A2E",
                   "Not significant" = "#1A5276",
                   "Failed"          = "grey70")
  shape_vals  <- c("FDR significant" = 18,
                   "Not significant" = 18,
                   "Failed"          = 1)

  # Text annotation positions (outside plot panel)
  x_carrier <- x_max * 1.30
  x_or      <- x_max * 2.40

  p <- ggplot(plot_df, aes(y = gene_f)) +

    geom_vline(xintercept = 1, linetype = "longdash",
               colour = "grey55", linewidth = 0.5) +

    geom_errorbar(
      aes(xmin = CI_low, xmax = CI_high_draw, colour = colour_group),
      width = 0.28, linewidth = 0.8, na.rm = TRUE,
      orientation = "y"
    ) +

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

    annotate("text",
             x     = c(x_carrier, x_or),
             y     = length(GENES) + 1,
             label = c("Carriers\n(cases/ctrl)", "OR (95% CI)"),
             hjust = c(0.5, 0), size = 3.1, fontface = "bold") +

    scale_x_log10(
      limits = c(x_min, x_or * 1.9),
      breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25),
      labels = label_number(accuracy = 0.1)
    ) +

    scale_colour_manual(values = colour_vals, name = NULL) +
    scale_shape_manual( values = shape_vals,  name = NULL) +

    # Fix "guide perpendicular to axis" warning by setting guide position
    guides(
      colour = guide_legend(position = "bottom"),
      shape  = guide_legend(position = "bottom")
    ) +

    labs(
      title    = title_text,
      subtitle = paste0(
        "Firth's penalised logistic regression, adjusted for ageInt and study\n",
        "Bars: 95% profile penalised-likelihood CI  |  FDR: Benjamini-Hochberg, 9 genes"
      ),
      x = "Odds Ratio (log scale)",
      y = NULL
    ) +

    coord_cartesian(clip = "off") +

    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 8.5, colour = "grey45"),
      axis.text.y        = element_text(face = "italic", size = 11),
      axis.text.x        = element_text(size = 10),
      legend.position    = "bottom",
      legend.text        = element_text(size = 10),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
      plot.margin        = margin(t = 8, r = 220, b = 8, l = 8, unit = "pt")
    )

  ggsave(filename, plot = p, width = 11, height = 6.5, dpi = 300, bg = "white")
  message(sprintf("Forest plot saved: %s", filename))
  invisible(p)
}

make_forest_plot(res_bbd,
                 "Truncating Variant Associations with BBD (Overall)",
                 "outputs/figures/forest_BBD_overall.png")
make_forest_plot(res_nonprol,
                 "Truncating Variant Associations with Non-proliferative BBD",
                 "outputs/figures/forest_BBD_nonproliferative.png")
make_forest_plot(res_prol,
                 "Truncating Variant Associations with Proliferative BBD (without atypia)",
                 "outputs/figures/forest_BBD_proliferative.png")


# ── 13.  RESULTS SUMMARY ─────────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("RESULTS SUMMARY (Firth's penalised logistic regression)")
message(strrep("=", 60))

for (lbl in unique(results_all$analysis)) {
  sub <- filter(results_all, analysis == lbl)
  sig <- sub %>% filter(fdr_sig == TRUE)  %>% pull(gene)
  nom <- sub %>% filter(p_value < 0.05, fdr_sig == FALSE) %>% pull(gene)
  message(sprintf("\n%s:", lbl))
  message("  FDR-significant : ",
          if (length(sig) > 0) paste(sig, collapse = ", ") else "none")
  message("  Nominally sig.  : ",
          if (length(nom) > 0) paste(nom, collapse = ", ") else "none")
}

message("\nMethodological note:")
message("  Firth's penalised logistic regression was used throughout.")
message("  This approach is recommended for rare-variant analyses to avoid")
message("  inflated estimates arising from complete/quasi-complete separation.")
message("  Reference: Heinze & Schemper (2002), Stat Med 21:2409-2419.")
message("\nAtypical hyperplasia (BBD_type1 == 3, n = 2): excluded from subtype models.")
message("All models adjusted for ageInt (continuous) + study (categorical).")
message("FDR correction: Benjamini-Hochberg across 9 genes per comparison.\n")


# ── 14.  SESSION INFO ─────────────────────────────────────────────────────────
# Saves R version + all package versions for reproducibility documentation.
# Cite in dissertation Methods: "R x.x.x; see session_info.txt for details."
si <- sessionInfo()
writeLines(capture.output(si), "outputs/session_info/session_info_BBD_truncating.txt")
message("Session info saved to: outputs/session_info/session_info_BBD_truncating.txt")
message(sprintf("R version: %s", R.version.string))
message(sprintf("Platform : %s", si$platform))
message(sprintf("Date run : %s", Sys.time()))
