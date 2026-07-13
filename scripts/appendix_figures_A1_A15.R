# =============================================================================
# Appendix A — Supplementary figures A1-A15
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil Population Health Sciences — 2025/26
#
# PURPOSE
#   Renders the fifteen Appendix A figures that support the Results chapter:
#     A1-A2   DCIS per-gene forest (truncating / missense)
#     A3-A4   BBD per-gene forest (n = 1,765 main genetic analysis)
#     A5-A6   LCIS per-gene forest
#     A7-A10  BBD subtype forests (non-proliferative n=431, proliferative n=274)
#     A11     Age distributions by cohort (QC)
#     A12     Top contributing studies by control count (QC)
#     A13     Ascertainment sensitivity of the DCIS truncating associations
#     A14     Leave-one-study-out analysis of the ATM-DCIS association
#     A15     BBD carrier depletion by subsequent cancer status (collider check)
#
#   Every panel is drawn from the SAVED, VERIFIED result tables (or, for the
#   QC panels A11/A12/A15, directly from the phenotype/genotype files) so a
#   figure can never disagree with a reported number. Appearance follows
#   scripts/_thesis_theme.R. Forest panels use the pre-revision (BRIDGES v17)
#   subtype coding, on which the per-gene subtype models are based.
#
# INPUTS : outputs/tables/{DCIS_LCIS_truncating_results, missense_results_FINAL,
#          BBD_truncating_FINAL, sensitivity_familial_population}.xlsx,
#          outputs/tables/loso_all_fits.rds, and the raw phenotype/genotype
#          files (A11/A12/A15 only; gitignored, present locally).
# OUTPUT : outputs/figures/appendix/figA1.png ... figA15.png (+ vector .pdf)
#
# DATA GOVERNANCE: outputs are non-disclosive summary statistics only.
# =============================================================================

suppressMessages({ library(ggplot2); library(dplyr); library(readxl); library(tidyr); library(readr); library(purrr) })
source("scripts/_thesis_theme.R")
OUT <- "outputs/figures/appendix"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")

# -- helpers ------------------------------------------------------------------
# Parse "OR (lo-hi)" text into numeric columns.
parse_ci <- function(x) {
  x <- gsub("[\u2013\u2212]", "-", as.character(x))
  m <- regmatches(x, regexec("([0-9.]+)\\s*\\(([0-9.]+)\\s*-\\s*([0-9.]+)\\)", x))
  t(sapply(m, function(v) if (length(v) == 4) as.numeric(v[2:4]) else c(NA, NA, NA)))
}

# One forest panel from a data frame with gene/OR/lo/hi (+ optional carriers).
forest_panel <- function(df, title, sig = character(0), file, xmax = NULL) {
  df <- df %>% mutate(gene = factor(gene, levels = rev(GENES)),
                      signif = gene %in% sig)
  p <- ggplot(df, aes(OR, gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = THESIS_PAL["rule"], linewidth = 0.4) +
    geom_errorbarh(aes(xmin = pmax(lo, 1e-2), xmax = hi, colour = signif), height = 0.25, linewidth = 0.6,
                   na.rm = TRUE) +
    geom_point(aes(colour = signif), size = 2.4, na.rm = TRUE) +
    { if (any(is.na(df$OR)))
        geom_text(data = df %>% filter(is.na(OR)), aes(x = 1, y = gene),
                  label = "not estimable", size = 2.5, fontface = "italic",
                  colour = THESIS_PAL["muted"]) } +
    scale_colour_manual(values = c(`TRUE` = unname(THESIS_PAL["signal"]),
                                   `FALSE` = "grey20"), guide = "none") +
    scale_x_log10() +
    labs(x = "Odds ratio (95% CI), log scale", y = NULL, title = title) +
    theme_thesis(grid = "x") +
    theme(axis.text.y = element_text(face = "italic"))
  if (!is.null(xmax)) p <- p + coord_cartesian(xlim = c(0.05, xmax))
  save_fig(p, file.path(OUT, file), width = 6.4, height = 4.0)
}

# -- A1/A2 DCIS, A5/A6 LCIS ---------------------------------------------------
read_forest <- function(path, sheet) {
  d <- read_excel(path, sheet)
  ci <- parse_ci(d[["OR (95% CI)"]])
  tibble(gene = d$Gene, OR = ci[,1], lo = ci[,2], hi = ci[,3],
         carr = suppressWarnings(as.integer(d[["Carriers (cases)"]])))
}
DL <- "outputs/tables/DCIS_LCIS_truncating_results.xlsx"
MS <- "outputs/tables/missense_results_FINAL.xlsx"
BB <- "outputs/tables/BBD_truncating_FINAL.xlsx"

forest_panel(read_forest(DL, "DCIS_truncating"),
  "Figure A1. DCIS, protein-truncating variants (n = 1,663)", c("ATM","BRCA2","CHEK2"), "figA1.png", 50)
forest_panel(read_forest(MS, "DCIS_missense"),
  "Figure A2. DCIS, missense variants (CADD >= 20; n = 1,663)", c("CHEK2","TP53"), "figA2.png", 12)
forest_panel(read_forest(BB, "BBD vs Controls"),
  "Figure A3. BBD, protein-truncating variants (n = 1,765)", character(0), "figA3.png", 2000)
forest_panel(read_forest(MS, "BBD_missense"),
  "Figure A4. BBD, missense variants (n = 1,765)", character(0), "figA4.png", 12)
forest_panel(read_forest(DL, "LCIS_truncating"),
  "Figure A5. LCIS, protein-truncating variants (n = 139)", character(0), "figA5.png", 300)
forest_panel(read_forest(MS, "LCIS_missense"),
  "Figure A6. LCIS, missense variants (n = 139)", character(0), "figA6.png", 50)
forest_panel(read_forest(BB, "Non-proliferative BBD"),
  "Figure A7. Non-proliferative BBD, truncating (n = 431)", character(0), "figA7.png", 500)
forest_panel(read_forest(BB, "Proliferative BBD"),
  "Figure A8. Proliferative BBD without atypia, truncating (n = 274)", character(0), "figA8.png", 300)
forest_panel(read_forest(MS, "BBD_nonprol_missense"),
  "Figure A9. Non-proliferative BBD, missense (n = 431)", character(0), "figA9.png", 50)
forest_panel(read_forest(MS, "BBD_prol_missense"),
  "Figure A10. Proliferative BBD without atypia, missense (n = 274)", character(0), "figA10.png", 50)

cat("A1-A10 forest panels written.\n")

# -- QC panels A11/A12/A15 need the raw files (gitignored, present locally) ----
PHENO <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC <- "concept_807_zhang_bridges_truncating.csv"
if (file.exists(PHENO)) {
  pheno <- suppressWarnings(read_delim(PHENO, delim = "\t", show_col_types = FALSE,
                                       na = c("", "NA", "888", "777", "999")))
  eur <- pheno %>% filter(ethnicityClass == 1, !is.na(ageInt))
  mc <- grep("Morphologygroup1", names(eur), value = TRUE)[1]

  # A11: age distributions
  cohorts <- bind_rows(
    eur %>% filter(status == 0, BBD_history == 1) %>% transmute(age = ageInt, grp = "BBD cases"),
    eur %>% filter(status == 0, BBD_history == 0) %>% transmute(age = ageInt, grp = "Cancer-free controls"),
    eur %>% filter(status == 2, .data[[mc]] == "Ductal") %>% transmute(age = ageInt, grp = "DCIS"),
    eur %>% filter(status == 2, .data[[mc]] == "Lobular") %>% transmute(age = ageInt, grp = "LCIS")) %>%
    mutate(grp = factor(grp, levels = c("BBD cases","Cancer-free controls","DCIS","LCIS")))
  means <- cohorts %>% group_by(grp) %>% summarise(m = mean(age), n = n(), .groups = "drop")
  labs11 <- setNames(sprintf("%s (n = %s; mean %.1f y)", means$grp, format(means$n, big.mark=","), means$m), means$grp)
  pA11 <- ggplot(cohorts, aes(age)) +
    geom_histogram(binwidth = 2.5, fill = THESIS_PAL["baseline"], colour = "white", linewidth = 0.2) +
    geom_vline(data = means, aes(xintercept = m), colour = THESIS_PAL["signal"], linewidth = 0.6, linetype = "dashed") +
    facet_wrap(~ grp, scales = "free_y", labeller = labeller(grp = labs11)) +
    labs(x = "Age at interview or diagnosis (years)", y = "Number of women",
         title = "Figure A11. Age distributions by cohort (quality-control check)") +
    theme_thesis(grid = "y")
  save_fig(pA11, file.path(OUT, "figA11.png"), width = 7.4, height = 5.2)

  # A12: top-20 studies by control count
  top <- eur %>% filter(status == 0) %>% count(study, sort = TRUE) %>% slice_head(n = 20)
  pA12 <- ggplot(top, aes(n, reorder(study, n))) +
    geom_col(fill = THESIS_PAL["baseline"]) +
    labs(x = "Number of cancer-free controls contributed", y = NULL,
         title = "Figure A12. Top 20 contributing studies (quality-control check)") +
    theme_thesis(grid = "x")
  save_fig(pA12, file.path(OUT, "figA12.png"), width = 6.4, height = 4.6)

  # A15: BBD carrier depletion by subsequent cancer status
  if (file.exists(TRUNC)) {
    trunc <- read_csv(TRUNC, show_col_types = FALSE, na = c("", "NA")) %>% select(-any_of("...1"))
    tcols <- paste0(GENES, "_truncating")
    carr <- trunc %>% transmute(BRIDGES_ID, across(all_of(tcols), ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))
    bh <- eur %>% filter(BBD_history == 1) %>% left_join(carr, by = "BRIDGES_ID") %>%
      mutate(across(all_of(tcols), ~ replace_na(.x, 0L)),
             grp = ifelse(status == 0, "Remained cancer-free", "Progressed to cancer"))
    dep <- bh %>% group_by(grp) %>%
      summarise(across(all_of(tcols), ~ mean(.x) * 100), .groups = "drop") %>%
      pivot_longer(-grp, names_to = "gene", values_to = "pct") %>%
      mutate(gene = factor(sub("_truncating", "", gene), levels = GENES),
             grp = factor(grp, levels = c("Remained cancer-free","Progressed to cancer")))
    pA15 <- ggplot(dep, aes(gene, pct, fill = grp)) +
      geom_col(position = "dodge", width = 0.75) +
      scale_fill_manual(values = c("Remained cancer-free" = unname(THESIS_PAL["baseline"]),
                                   "Progressed to cancer" = unname(THESIS_PAL["signal"])), name = NULL) +
      labs(x = NULL, y = "Truncating-variant carrier frequency (%)",
           title = "Figure A15. BBD carrier depletion by subsequent cancer status") +
      theme_thesis(grid = "y") + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1))
    save_fig(pA15, file.path(OUT, "figA15.png"), width = 6.8, height = 4.2)
  }
  cat("A11/A12/A15 QC panels written.\n")
} else {
  message("Raw phenotype file not present; skipping A11/A12/A15 (QC panels).")
}

# -- A13 ascertainment sensitivity --------------------------------------------
fam <- read_excel("outputs/tables/sensitivity_familial_population.xlsx")
spec_cols <- c("(A) All studies [PRIMARY]", "(B) Population only", "(C) famHist-adjusted")
a13 <- fam %>% filter(Variant == "truncating", Gene %in% c("ATM","CHEK2","BRCA2")) %>%
  select(Gene, all_of(spec_cols)) %>%
  pivot_longer(-Gene, names_to = "spec", values_to = "cell")
ci <- parse_ci(a13$cell)
a13 <- a13 %>% mutate(OR = ci[,1], lo = ci[,2], hi = ci[,3],
                      spec = factor(spec, levels = spec_cols,
                                    labels = c("All studies","Population only","Family-history adjusted")),
                      gene = factor(Gene, levels = c("BRCA2","CHEK2","ATM"))) %>%
  filter(!is.na(OR))
pA13 <- ggplot(a13, aes(OR, gene, colour = spec)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = THESIS_PAL["rule"], linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = THESIS_GROUPS, name = NULL) +
  scale_x_log10() +
  labs(x = "Odds ratio (95% CI), log scale", y = NULL,
       title = "Figure A13. Ascertainment sensitivity of the DCIS truncating associations") +
  theme_thesis(grid = "x") + theme(axis.text.y = element_text(face = "italic"))
save_fig(pA13, file.path(OUT, "figA13.png"), width = 6.6, height = 4.0)

# -- A14 leave-one-study-out (ATM-DCIS truncating) ----------------------------
lo <- readRDS("outputs/tables/loso_all_fits.rds")
full_or <- lo %>% filter(gene == "ATM", vtype == "truncating", grepl("full data", dropped)) %>% pull(OR)
a14 <- lo %>% filter(gene == "ATM", vtype == "truncating", !grepl("full data", dropped), !is.na(OR)) %>%
  arrange(OR) %>% mutate(dropped = factor(dropped, levels = dropped))
extremes <- c(as.character(a14$dropped[1]), as.character(a14$dropped[nrow(a14)]))
pA14 <- ggplot(a14, aes(OR, dropped)) +
  annotate("rect", xmin = min(a14$OR), xmax = max(a14$OR), ymin = -Inf, ymax = Inf,
           fill = THESIS_PAL["muted"], alpha = 0.25) +
  geom_vline(xintercept = full_or, colour = THESIS_PAL["signal"], linewidth = 0.9) +
  geom_point(aes(colour = dropped %in% extremes), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = unname(THESIS_PAL["signal"]), `FALSE` = "grey20"), guide = "none") +
  labs(x = "ATM-DCIS odds ratio when that study is removed", y = NULL,
       title = "Figure A14. Leave-one-study-out analysis of the ATM-DCIS association",
       caption = sprintf("%d leave-one-out refits; OR range %.2f-%.2f; full-data OR %.2f (red line); all p < 0.05.",
                         nrow(a14), min(a14$OR), max(a14$OR), full_or)) +
  theme_thesis(grid = "x")
save_fig(pA14, file.path(OUT, "figA14.png"), width = 6.8, height = 5.0)

cat("A13/A14 written. Appendix A complete: 15 figures in", OUT, "\n")
