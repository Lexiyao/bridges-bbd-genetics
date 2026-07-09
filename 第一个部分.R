# =============================================================================
# Dissertation Analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil — Objective 1: BBD Truncating Variant Analysis
#
# COLUMN NAMES (verified from actual BRIDGES dataset, not dictionary):
#   BRIDGES_ID   — participant ID (join key)
#   status       — 0=control, 1=invasive, 2=in-situ, 3=unknown invasiveness
#   ageInt       — age at interview (continuous)
#   study        — BCAC study acronym (categorical)
#   BBD_history  — 1=yes, 0=no (within status==0 only)
#   BBD_type1    — 1=non-proliferative, 2=proliferative w/o atypia,
#                  3=atypical hyperplasia (n=2, excluded)
#
# CORRECT BBD DEFINITION:
#   Cases    = status == 0 AND BBD_history == 1  (expected n = 1,765)
#   Controls = status == 0 AND BBD_history == 0  (expected n = 12,715)
#   status == 3 ("case of unknown invasiveness") is NOT BBD — excluded
#
# Nine genes: BRCA1, BRCA2, PALB2, CHEK2, ATM, BARD1, RAD51C, RAD51D, TP53
# Model: outcome ~ gene_carrier + ageInt + study
# FDR: Benjamini-Hochberg across 9 genes per comparison
# =============================================================================


# ── 0.  PACKAGES ──────────────────────────────────────────────────────────────
required <- c("tidyverse", "broom", "writexl", "ggplot2", "scales")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(tidyverse)
library(broom)
library(writexl)
library(ggplot2)
library(scales)


# ── 1.  FILE PATHS ────────────────────────────────────────────────────────────
PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
OUT_XLSX   <- "BBD_truncating_results.xlsx"


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
  TRUNC_FILE, show_col_types = FALSE,
  na = c("", "NA")
) %>%
  # The CSV exports an unnamed row-index column read as "...1" — drop it
  select(-any_of("...1"))

message(sprintf("Phenotype rows : %d", nrow(pheno_raw)))
message(sprintf("Truncating rows: %d", nrow(trunc_raw)))


# ── 4.  VALIDATE COLUMNS ─────────────────────────────────────────────────────
# Using the actual column names observed in the BRIDGES dataset
required_pheno <- c("BRIDGES_ID", "status", "ageInt", "study",
                    "BBD_history", "BBD_type1")

missing_pheno <- setdiff(required_pheno, names(pheno_raw))
if (length(missing_pheno) > 0) {
  stop("Missing phenotype columns: ", paste(missing_pheno, collapse = ", "),
       "\nAvailable: ", paste(names(pheno_raw), collapse = ", "))
}
message("All required phenotype columns present [OK]")

genes_found   <- intersect(GENES, names(trunc_raw))
genes_missing <- setdiff(GENES, names(trunc_raw))
if (length(genes_missing) > 0)
  warning("Gene columns not found in truncating file: ",
          paste(genes_missing, collapse = ", "))
message("Gene columns found: ", paste(genes_found, collapse = ", "))

if (!"BRIDGES_ID" %in% names(trunc_raw))
  stop("BRIDGES_ID not found in truncating file. ",
       "Available: ", paste(names(trunc_raw), collapse = ", "))


# ── 5.  MERGE ─────────────────────────────────────────────────────────────────
dat_full <- pheno_raw %>%
  left_join(
    trunc_raw %>% select(BRIDGES_ID, all_of(genes_found)),
    by = "BRIDGES_ID"
  ) %>%
  # Recode NA in carrier columns to 0 (non-carrier) — standard BRIDGES practice
  mutate(across(all_of(genes_found), ~ replace_na(.x, 0L)))

message(sprintf("Merged dataset: %d rows", nrow(dat_full)))


# ── 6.  DEFINE BBD COHORT ─────────────────────────────────────────────────────
# Restrict to status == 0 first, then split by BBD_history.
# status == 3 is explicitly excluded (unknown invasiveness per data dictionary).

dat_bbd_all <- dat_full %>%
  filter(status == 0, !is.na(BBD_history)) %>%
  mutate(outcome = as.integer(BBD_history == 1))

n_bbd  <- sum(dat_bbd_all$outcome == 1)
n_ctrl <- sum(dat_bbd_all$outcome == 0)

message(sprintf("\nBBD cases   (status==0, BBD_history==1): n = %d  [expected 1,765]",  n_bbd))
message(sprintf("Controls    (status==0, BBD_history==0): n = %d  [expected 12,715]", n_ctrl))

if (abs(n_bbd  - 1765)  > 100) warning("BBD case count deviates from expected 1,765.")
if (abs(n_ctrl - 12715) > 200) warning("Control count deviates from expected 12,715.")


# ── 7.  BBD SUBTYPES ──────────────────────────────────────────────────────────
# Non-proliferative (BBD_type1 == 1) vs all controls
dat_nonprol <- dat_bbd_all %>%
  filter(outcome == 0 | (outcome == 1 & BBD_type1 == 1)) %>%
  mutate(outcome = as.integer(BBD_type1 == 1))

# Proliferative without atypia (BBD_type1 == 2) vs all controls
dat_prol <- dat_bbd_all %>%
  filter(outcome == 0 | (outcome == 1 & BBD_type1 == 2)) %>%
  mutate(outcome = as.integer(BBD_type1 == 2))

# Atypical hyperplasia (BBD_type1 == 3): n = 2, excluded from all models

message(sprintf("\nNon-proliferative cases:        n = %d  [expected ~431]",
                sum(dat_nonprol$outcome == 1)))
message(sprintf("Proliferative w/o atypia cases: n = %d  [expected ~274]",
                sum(dat_prol$outcome == 1)))
message("Atypical hyperplasia (n = 2): excluded from subtype models")


# ── 8.  LOGISTIC REGRESSION HELPER ────────────────────────────────────────────
# Model: logit(P(outcome)) = β0 + β_carrier + β_age*ageInt + Σ β_j*study_j
# - ageInt: continuous (years)
# - study: unordered factor (R creates K-1 dummy variables automatically)
# - Profile-likelihood 95% CIs used (more reliable than Wald for rare variants)

run_logistic <- function(data, gene) {
  
  n_carrier_case <- sum(data[[gene]] == 1 & data$outcome == 1, na.rm = TRUE)
  n_carrier_ctrl <- sum(data[[gene]] == 1 & data$outcome == 0, na.rm = TRUE)
  n_case         <- sum(data$outcome == 1, na.rm = TRUE)
  n_ctrl         <- sum(data$outcome == 0, na.rm = TRUE)
  
  if ((n_carrier_case + n_carrier_ctrl) == 0) {
    message(sprintf("  Skipping %s — no carriers observed in this comparison", gene))
    return(tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, converged = NA))
  }
  
  model_data <- data %>%
    select(outcome, carrier = all_of(gene), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()
  
  fit <- tryCatch(
    glm(outcome ~ carrier + ageInt + study,
        data = model_data, family = binomial(link = "logit")),
    error = function(e) {
      message(sprintf("  Model failed for %s: %s", gene, conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(fit)) {
    return(tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }
  
  coef_row <- tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "carrier")
  
  tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
         OR        = coef_row$estimate,
         CI_low    = coef_row$conf.low,
         CI_high   = coef_row$conf.high,
         p_value   = coef_row$p.value,
         converged = fit$converged)
}


# ── 9.  RUN ALL ANALYSES + FDR ────────────────────────────────────────────────
run_panel <- function(data, label) {
  message(sprintf("\nRunning: %s", label))
  map_dfr(GENES, ~ run_logistic(data, .x)) %>%
    mutate(
      analysis = label,
      p_fdr    = p.adjust(p_value, method = "BH"),
      fdr_sig  = p_fdr < 0.05
    ) %>%
    select(analysis, gene, n_case, n_ctrl,
           n_carrier_case, n_carrier_ctrl,
           OR, CI_low, CI_high, p_value, p_fdr, fdr_sig, converged)
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


# ── 12.  FOREST PLOT FUNCTION ─────────────────────────────────────────────────
make_forest_plot <- function(res_df, title_text, filename) {
  
  plot_df <- res_df %>%
    mutate(
      gene_f        = factor(gene, levels = rev(GENES)),
      colour_group  = case_when(
        is.na(OR) ~ "No carriers",
        fdr_sig   ~ "FDR significant",
        TRUE      ~ "Not significant"
      ),
      or_label      = if_else(
        is.na(OR), "No carriers",
        sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)
      ),
      carrier_label = sprintf("%d / %d", n_carrier_case, n_carrier_ctrl)
    )
  
  x_min <- min(c(plot_df$CI_low,  0.5), na.rm = TRUE) * 0.75
  x_max <- max(c(plot_df$CI_high, 3.0), na.rm = TRUE) * 1.40
  
  colour_vals <- c("FDR significant" = "#B03A2E",
                   "Not significant" = "#1A5276",
                   "No carriers"     = "grey70")
  shape_vals  <- c("FDR significant" = 18,
                   "Not significant" = 18,
                   "No carriers"     = 1)
  
  p <- ggplot(plot_df, aes(y = gene_f)) +
    
    geom_vline(xintercept = 1, linetype = "longdash",
               colour = "grey55", linewidth = 0.5) +
    
    geom_errorbarh(
      aes(xmin = CI_low, xmax = CI_high, colour = colour_group),
      height = 0.28, linewidth = 0.8, na.rm = TRUE
    ) +
    
    geom_point(
      aes(x = OR, colour = colour_group, shape = colour_group),
      size = 4, na.rm = TRUE
    ) +
    
    geom_text(
      aes(x = x_max * 1.10, label = carrier_label),
      hjust = 0.5, size = 3, colour = "grey25"
    ) +
    
    geom_text(
      aes(x = x_max * 1.55, label = or_label),
      hjust = 0, size = 3, colour = "grey15"
    ) +
    
    annotate("text",
             x     = c(x_max * 1.10, x_max * 1.55),
             y     = length(GENES) + 1,
             label = c("Carriers\n(cases/ctrl)", "OR (95% CI)"),
             hjust = c(0.5, 0), size = 3.1, fontface = "bold") +
    
    scale_x_log10(
      limits = c(x_min, x_max * 2.8),
      breaks = c(0.25, 0.5, 1, 2, 5, 10, 25),
      labels = label_number(accuracy = 0.1)
    ) +
    
    scale_colour_manual(values = colour_vals, name = NULL) +
    scale_shape_manual(values  = shape_vals,  name = NULL) +
    
    labs(
      title    = title_text,
      subtitle = paste0(
        "Binary logistic regression, adjusted for ageInt (continuous) and study (categorical)\n",
        "Bars: 95% profile-likelihood CI  |  FDR correction: Benjamini-Hochberg across 9 genes"
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
      plot.margin        = margin(t = 8, r = 200, b = 8, l = 8, unit = "pt")
    )
  
  ggsave(filename, plot = p, width = 11, height = 6.5, dpi = 300, bg = "white")
  message(sprintf("Forest plot saved: %s", filename))
  invisible(p)
}

make_forest_plot(res_bbd,
                 "Truncating Variant Associations with BBD (Overall)",
                 "forest_BBD_overall.png")

make_forest_plot(res_nonprol,
                 "Truncating Variant Associations with Non-proliferative BBD",
                 "forest_BBD_nonproliferative.png")

make_forest_plot(res_prol,
                 "Truncating Variant Associations with Proliferative BBD (without atypia)",
                 "forest_BBD_proliferative.png")


# ── 13.  RESULTS SUMMARY ─────────────────────────────────────────────────────
message("\n", strrep("=", 60))
message("RESULTS SUMMARY")
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

message("\nNote: Atypical hyperplasia (BBD_type1 == 3, n = 2) excluded from ",
        "all subtype models — insufficient sample size for stable estimates.")
message("All models adjusted for ageInt (continuous) + study (categorical).")
message("FDR correction applied within each comparison across 9 genes.\n")