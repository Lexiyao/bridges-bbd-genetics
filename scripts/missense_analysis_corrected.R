# =============================================================================
# Dissertation Analysis — Missense Variant Analysis (CORRECTED FINAL VERSION)
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil
#
# DATA GOVERNANCE
#   Data accessed under BCAC/BRIDGES approved secondary use agreement.
#   All analysis conducted on non-identifiable summary-level outputs only.
#   Raw data files are not redistributed or committed to version control.
#
# GENOTYPE QC PROVENANCE
#   Variant-level QC applied upstream by BRIDGES consortium prior to release.
#   Reference: Dorling et al. (2021) NEJM 384:428-439 (BRIDGES QC protocol).
#   Missense carriers defined by GENE_CADD.phred.01 == 1, a pre-computed
#   BRIDGES binary indicator for CADD phred score >= 20 (top 1% of all human
#   genome SNPs by predicted deleteriousness; Kircher et al. 2014 Nat Genet).
#   Despite the "01" suffix suggesting a quintile, verified from data that the
#   threshold is a fixed CADD phred >= 20 cutoff, not a quintile boundary.
#   No additional variant filtering applied here.
#   Principal components (PCs) not available in phenotype file v17.
#   Ancestry controlled by ethnicityClass == 1 + study adjustment.
#   This is a documented limitation — see dissertation Section 5.4.
#
# COMPUTATIONAL ENVIRONMENT
#   R version and package versions saved to session_info_missense.txt at end.
#
# KEY CORRECTIONS vs previous version:
#   1. ethnicityClass == 1 filter applied to ALL cohorts (European ancestry)
#      This brings BBD n: 2019 → 1765, controls: 52227 → 34424
#   2. DCIS/LCIS use MorphologygroupIndex_corr (text labels), not ICD-O codes
#      DCIS = "Ductal" (n=1663), LCIS = "Lobular" (n=139)
#   3. Forest plot text overlap fixed (wider spacing)
#
# CARRIER DEFINITION: GENE_CADD.phred.01 == 1 (CADD phred >= 20, top 1% most
#   deleterious variants genome-wide; verified from data — NOT a quintile cutoff)
# MODEL: Firth's penalised logistic regression
# COVARIATES: ageInt (continuous) + study (categorical)
# FDR: Benjamini-Hochberg across 9 genes per comparison
# =============================================================================

# ── 0.  PACKAGES ──────────────────────────────────────────────────────────────
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


# ── 1.  FILE PATHS ────────────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
MISS_FILE  <- "concept_807_zhang_bridges_missense.csv"
OUT_XLSX   <- "outputs/tables/missense_results_FINAL.xlsx"

# Ensure output directories exist (script previously wrote PNGs/xlsx to the
# project root; now writes to the organised outputs/ tree).
dir.create("outputs/tables",       recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures",      recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/session_info", recursive = TRUE, showWarnings = FALSE)

# Working-directory guard: scripts use paths relative to the project root.
for (f in c(PHENO_FILE, MISS_FILE)) {
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

message("Loading missense file...")
miss_raw <- read_csv(MISS_FILE, show_col_types = FALSE,
                     na = c("", "NA")) %>%
  select(-any_of("...1"))


# ── 4.  AUTO-DETECT CADD.phred.01 COLUMNS ────────────────────────────────────
detect_cadd_col <- function(gene, col_names) {
  m <- col_names[grepl(paste0("^", gene, "_CADD\\.phred\\.01$"),
                        col_names, ignore.case = TRUE)]
  if (length(m) > 0) return(m[1])
  m2 <- col_names[grepl(paste0("^", gene, ".*CADD.*01"), col_names,
                          ignore.case = TRUE)]
  if (length(m2) > 0) return(m2[1])
  NA_character_
}

cadd_col_map <- setNames(
  sapply(GENES, detect_cadd_col, col_names = names(miss_raw)), GENES)
genes_found <- names(cadd_col_map)[!is.na(cadd_col_map)]

message("CADD columns mapped: ", paste(genes_found, collapse = ", "))


# ── 5.  MERGE ─────────────────────────────────────────────────────────────────
miss_select <- miss_raw %>%
  select(BRIDGES_ID, all_of(unname(cadd_col_map[genes_found]))) %>%
  rename(!!!setNames(cadd_col_map[genes_found], genes_found))

# Join-coverage guard: every phenotype participant must have a genotype record,
# else left_join + replace_na(0) would silently code them as non-carriers.
n_unmatched <- pheno_raw %>%
  anti_join(miss_select, by = "BRIDGES_ID") %>%
  nrow()
if (n_unmatched > 0)
  stop(sprintf("%d phenotype IDs have no missense-genotype record and would be ",
               n_unmatched),
       "silently coded as non-carriers. Aborting.", call. = FALSE)

dat_full <- pheno_raw %>%
  left_join(miss_select, by = "BRIDGES_ID") %>%
  mutate(across(all_of(genes_found), ~ replace_na(as.integer(.x), 0L)))

message(sprintf("Merged dataset: %d rows", nrow(dat_full)))


# ── 6.  DEFINE COHORTS (ALL with ethnicityClass == 1) ─────────────────────────
#
# CRITICAL: all cohorts restricted to European ancestry (ethnicityClass == 1)
# This is consistent with the previously completed DCIS/LCIS truncating analyses
# and resolves the BBD count discrepancy (2019 → 1765).

# 6a. BBD overall
dat_bbd_all <- dat_full %>%
  filter(status == 0,
         ethnicityClass == 1,        # European ancestry
         !is.na(ageInt),             # age covariate present (matches truncating/QC cohort)
         !is.na(BBD_history)) %>%
  mutate(outcome = as.integer(BBD_history == 1))

# 6b. BBD subtypes
dat_nonprol <- dat_bbd_all %>%
  filter(outcome == 0 | (outcome == 1 & BBD_type1 == 1)) %>%
  mutate(outcome = if_else(BBD_type1 == 1, 1L, 0L, missing = 0L))

dat_prol <- dat_bbd_all %>%
  filter(outcome == 0 | (outcome == 1 & BBD_type1 == 2)) %>%
  mutate(outcome = if_else(BBD_type1 == 2, 1L, 0L, missing = 0L))

# 6c. DCIS and LCIS
# Use MorphologygroupIndex_corr (text labels) — consistent with dat_clean
# 888 = missing morphology → excluded
# Mixed ductal/lobular → excluded (n=63), consistent with previous analysis

# !is.na(ageInt) added to every cohort so reported N and carrier counts match
# the model cohort (drop_na() in run_firth removes missing-age rows). Without
# it, reported controls were 37,568 (vs the correct 34,424) while the model
# silently used 34,424 — an internal inconsistency.
controls_pool <- dat_full %>%
  filter(status == 0, ethnicityClass == 1, !is.na(ageInt)) %>%
  mutate(outcome = 0L)

dcis_cases <- dat_full %>%
  filter(status == 2,
         ethnicityClass == 1,
         !is.na(ageInt),
         MorphologygroupIndex_corr == "Ductal") %>%
  mutate(outcome = 1L)

lcis_cases <- dat_full %>%
  filter(status == 2,
         ethnicityClass == 1,
         !is.na(ageInt),
         MorphologygroupIndex_corr == "Lobular") %>%
  mutate(outcome = 1L)

dat_dcis <- bind_rows(controls_pool, dcis_cases)
dat_lcis <- bind_rows(controls_pool, lcis_cases)

# Report cohort sizes — verify against previous analyses
message("\n── Cohort sizes (ethnicityClass == 1 applied) ──")
message(sprintf("  BBD overall          : %d cases / %d controls  [expected 1765 / 12715]",
                sum(dat_bbd_all$outcome == 1), sum(dat_bbd_all$outcome == 0)))
message(sprintf("  Non-proliferative BBD: %d cases / %d controls  [expected 431]",
                sum(dat_nonprol$outcome == 1), sum(dat_nonprol$outcome == 0)))
message(sprintf("  Proliferative BBD    : %d cases / %d controls  [expected 274]",
                sum(dat_prol$outcome == 1),    sum(dat_prol$outcome == 0)))
message(sprintf("  DCIS                 : %d cases / %d controls  [expected 1663 / 34424]",
                sum(dat_dcis$outcome == 1),    sum(dat_dcis$outcome == 0)))
message(sprintf("  LCIS                 : %d cases / %d controls  [expected 139 / 34424]",
                sum(dat_lcis$outcome == 1),    sum(dat_lcis$outcome == 0)))


# ── 7.  FIRTH'S LOGISTIC REGRESSION HELPER ────────────────────────────────────
run_firth <- function(data, gene) {
  n_carrier_case <- sum(data[[gene]] == 1 & data$outcome == 1, na.rm = TRUE)
  n_carrier_ctrl <- sum(data[[gene]] == 1 & data$outcome == 0, na.rm = TRUE)
  n_case         <- sum(data$outcome == 1, na.rm = TRUE)
  n_ctrl         <- sum(data$outcome == 0, na.rm = TRUE)

  if ((n_carrier_case + n_carrier_ctrl) == 0) {
    message(sprintf("  %s: 0 carriers — skipped", gene))
    return(tibble(gene, n_case, n_ctrl,
                  n_carrier_case = 0L, n_carrier_ctrl = 0L,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "no carriers"))
  }

  model_data <- data %>%
    select(outcome, carrier = all_of(gene), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()

  fit <- tryCatch(
    logistf(outcome ~ carrier + ageInt + study,
            data = model_data, firth = TRUE, pl = TRUE),
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

  idx <- which(names(coef(fit)) == "carrier")
  tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
         OR      = exp(coef(fit)[idx]),
         CI_low  = exp(fit$ci.lower[idx]),
         CI_high = exp(fit$ci.upper[idx]),
         p_value = fit$prob[idx],
         method  = "Firth")
}


# ── 8.  RUN ALL ANALYSES + FDR ────────────────────────────────────────────────
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

res_miss_bbd     <- run_panel(dat_bbd_all, "BBD vs Controls (missense)")
res_miss_nonprol <- run_panel(dat_nonprol, "Non-proliferative BBD vs Controls (missense)")
res_miss_prol    <- run_panel(dat_prol,    "Proliferative BBD vs Controls (missense)")
res_miss_dcis    <- run_panel(dat_dcis,    "DCIS vs Controls (missense)")
res_miss_lcis    <- run_panel(dat_lcis,    "LCIS vs Controls (missense)")

results_miss_all <- bind_rows(res_miss_bbd, res_miss_nonprol, res_miss_prol,
                               res_miss_dcis, res_miss_lcis)


# ── 9.  FORMAT TABLE ─────────────────────────────────────────────────────────
fmt_table <- function(df) {
  df %>%
    mutate(
      `OR (95% CI)` = case_when(
        is.na(OR) ~ "—",
        TRUE      ~ sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)),
      `p-value` = case_when(
        is.na(p_value)  ~ "—",
        p_value < 0.001 ~ "<0.001",
        TRUE            ~ sprintf("%.3f", p_value)),
      `p-FDR (BH)` = case_when(
        is.na(p_fdr)  ~ "—",
        p_fdr < 0.001 ~ "<0.001",
        TRUE          ~ sprintf("%.3f", p_fdr)),
      `FDR sig.` = if_else(fdr_sig, "Yes*", "No", missing = "—")
    ) %>%
    select(Gene = gene, `N (cases)` = n_case, `N (controls)` = n_ctrl,
           `Carriers (cases)` = n_carrier_case,
           `Carriers (controls)` = n_carrier_ctrl,
           `OR (95% CI)`, `p-value`, `p-FDR (BH)`, `FDR sig.`)
}

for (lbl in unique(results_miss_all$analysis)) {
  cat(strrep("─", 70), "\n", lbl, "\n", strrep("─", 70), "\n")
  print(fmt_table(filter(results_miss_all, analysis == lbl)), n = Inf)
  cat("\n")
}


# ── 10.  EXPORT ──────────────────────────────────────────────────────────────
write_xlsx(
  list("BBD_missense"         = fmt_table(res_miss_bbd),
       "BBD_nonprol_missense" = fmt_table(res_miss_nonprol),
       "BBD_prol_missense"    = fmt_table(res_miss_prol),
       "DCIS_missense"        = fmt_table(res_miss_dcis),
       "LCIS_missense"        = fmt_table(res_miss_lcis),
       "All_raw"              = results_miss_all),
  path = OUT_XLSX)
message(sprintf("\nExported to: %s", OUT_XLSX))


# ── 11.  FOREST PLOT (text overlap fixed) ────────────────────────────────────
# NOTE: the canonical dissertation figures are now regenerated centrally — in one
# unified house style — by scripts/build_exhibits.R via the shared forest_plot()
# helper (scripts/_exhibits_helpers.R), drawn from the saved result tables. This
# inline version is retained only for standalone runs of this script.
make_forest_plot <- function(res_df, title_text, filename) {

  plot_df <- res_df %>%
    mutate(
      gene_f        = factor(gene, levels = rev(GENES)),
      colour_group  = case_when(
        is.na(OR) ~ "No carriers",
        fdr_sig   ~ "FDR significant",
        TRUE      ~ "Not significant"),
      # Truncate very wide CIs for display only (does not affect statistics)
      CI_high_plot  = pmin(CI_high, 50),
      or_label      = if_else(is.na(OR), "No carriers",
                               sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)),
      carrier_label = sprintf("%d / %d", n_carrier_case, n_carrier_ctrl)
    )

  # Fixed text columns — not dependent on CI width (fixes overlap bug)
  x_plot_max <- 50    # display limit
  x_carrier  <- 80    # fixed position for carrier count column
  x_or       <- 120   # fixed position for OR (95% CI) column

  colour_vals <- c("FDR significant" = "#B03A2E",
                   "Not significant" = "#1A5276",
                   "No carriers"     = "grey70")
  shape_vals  <- c("FDR significant" = 18,
                   "Not significant" = 18,
                   "No carriers"     = 1)

  p <- ggplot(plot_df, aes(y = gene_f)) +

    geom_vline(xintercept = 1, linetype = "longdash",
               colour = "grey55", linewidth = 0.5) +

    # CI bars capped at x_plot_max for readability
    geom_errorbar(
      aes(xmin = pmax(CI_low, 0.05), xmax = CI_high_plot,
          colour = colour_group),
      width = 0.28, linewidth = 0.8, na.rm = TRUE, orientation = "y"
    ) +

    geom_point(aes(x = OR, colour = colour_group, shape = colour_group),
               size = 4, na.rm = TRUE) +

    # Separator line between plot and text
    geom_vline(xintercept = 65, colour = "grey85", linewidth = 0.4) +

    # Carrier count column
    geom_text(aes(x = x_carrier, label = carrier_label),
              hjust = 0.5, size = 3, colour = "grey25") +

    # OR (95% CI) column
    geom_text(aes(x = x_or, label = or_label),
              hjust = 0, size = 3, colour = "grey15") +

    # Column headers
    annotate("text",
             x     = c(x_carrier, x_or),
             y     = length(GENES) + 1,
             label = c("Carriers\n(cases/ctrl)", "OR (95% CI)"),
             hjust = c(0.5, 0), size = 3.1, fontface = "bold") +

    scale_x_log10(
      limits = c(0.05, 300),
      breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25, 50),
      labels = label_number(accuracy = 0.1)
    ) +

    scale_colour_manual(values = colour_vals, name = NULL) +
    scale_shape_manual( values = shape_vals,  name = NULL) +
    guides(colour = guide_legend(position = "bottom"),
           shape  = guide_legend(position = "bottom")) +

    labs(
      title    = title_text,
      subtitle = paste0(
        "Firth's penalised logistic regression | European ancestry (ethnicityClass == 1)\n",
        "Adjusted for ageInt (continuous) and study (categorical)\n",
        "Carriers: CADD phred ≥20 (top 1% deleterious, GENE_CADD.phred.01==1)  |  FDR: BH, 9 genes"
      ),
      x = "Odds Ratio (log scale)", y = NULL
    ) +

    coord_cartesian(clip = "off") +

    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 7.5, colour = "grey45"),
      axis.text.y        = element_text(face = "italic", size = 11),
      axis.text.x        = element_text(size = 10),
      legend.position    = "bottom",
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
      plot.margin        = margin(t = 8, r = 260, b = 8, l = 8, unit = "pt")
    )

  ggsave(filename, plot = p, width = 13, height = 6.5, dpi = 300, bg = "white")
  message(sprintf("Saved: %s", filename))
  invisible(p)
}

make_forest_plot(res_miss_bbd,
                 "Missense Variant Associations with BBD (Overall)",
                 "outputs/figures/forest_missense_BBD_corrected.png")
make_forest_plot(res_miss_nonprol,
                 "Missense Variant Associations with Non-proliferative BBD",
                 "outputs/figures/forest_missense_BBD_nonprol_corrected.png")
make_forest_plot(res_miss_prol,
                 "Missense Variant Associations with Proliferative BBD (without atypia)",
                 "outputs/figures/forest_missense_BBD_prol_corrected.png")
make_forest_plot(res_miss_dcis,
                 "Missense Variant Associations with DCIS",
                 "outputs/figures/forest_missense_DCIS_corrected.png")
make_forest_plot(res_miss_lcis,
                 "Missense Variant Associations with LCIS",
                 "outputs/figures/forest_missense_LCIS_corrected.png")


# ── 12.  RESULTS SUMMARY ─────────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("MISSENSE RESULTS SUMMARY (corrected, ethnicityClass==1)")
message(strrep("=", 60))

for (lbl in unique(results_miss_all$analysis)) {
  sub <- filter(results_miss_all, analysis == lbl)
  sig <- sub %>% filter(fdr_sig == TRUE)  %>% pull(gene)
  nom <- sub %>% filter(p_value < 0.05, fdr_sig == FALSE) %>% pull(gene)
  message(sprintf("\n%s:", lbl))
  message("  FDR-sig  : ", if(length(sig)>0) paste(sig,collapse=", ") else "none")
  message("  Nom. sig : ", if(length(nom)>0) paste(nom,collapse=", ") else "none")
}

message("\nNote: ethnicityClass == 1 applied to all cohorts.")
message("DCIS defined as MorphologygroupIndex_corr == 'Ductal'")
message("LCIS defined as MorphologygroupIndex_corr == 'Lobular'")
message("888 (missing morphology), Mixed, Other, Papillary excluded.\n")


# ── 13.  SESSION INFO ─────────────────────────────────────────────────────────
si <- sessionInfo()
writeLines(capture.output(si), "outputs/session_info/session_info_missense.txt")
message("Session info saved to: outputs/session_info/session_info_missense.txt")
message(sprintf("R version: %s", R.version.string))
message(sprintf("Platform : %s", si$platform))
message(sprintf("Date run : %s", Sys.time()))
