# =============================================================================
# Build all gene-panel forest figures from the SAVED result tables
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# Reads the verified result workbooks (outputs/tables/*.xlsx) and redraws every
# forest plot through the single shared forest_plot() helper, so all figures
# share one style and every value matches the analysis output exactly. No model
# is re-fitted here; this script only visualises numbers that already exist.
#
# Run:  Rscript scripts/build_exhibits.R
# =============================================================================

suppressPackageStartupMessages({library(readxl); library(dplyr)})
source("scripts/_exhibits_helpers.R")
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# One row per figure: which workbook/sheet, which `analysis` slice, the output
# filename (matching existing names so document references stay valid), title,
# variant class (selects the caption) and the x-axis cap.
SPEC <- tribble(
  ~file,                                          ~sheet,         ~analysis,                                              ~out,                                       ~title,                                                       ~class,        ~cap,
  "outputs/tables/DCIS_LCIS_truncating_results.xlsx","Raw_results","DCIS vs Controls (truncating)",                        "forest_DCIS_truncating.png",               "Truncating Variant Associations with DCIS",                  "truncating",   25,
  "outputs/tables/DCIS_LCIS_truncating_results.xlsx","Raw_results","LCIS vs Controls (truncating)",                        "forest_LCIS_truncating.png",               "Truncating Variant Associations with LCIS",                  "truncating",   25,
  "outputs/tables/BBD_truncating_FINAL.xlsx",        "Raw results","BBD vs Controls",                                      "forest_BBD_overall.png",                   "Truncating Variant Associations with BBD (Overall)",         "truncating",   25,
  "outputs/tables/BBD_truncating_FINAL.xlsx",        "Raw results","Non-proliferative BBD vs Controls",                    "forest_BBD_nonproliferative.png",          "Truncating Variant Associations with Non-proliferative BBD", "truncating",   25,
  "outputs/tables/BBD_truncating_FINAL.xlsx",        "Raw results","Proliferative BBD vs Controls (without atypia)",       "forest_BBD_proliferative.png",             "Truncating Variant Associations with Proliferative BBD",     "truncating",   25,
  "outputs/tables/missense_results_FINAL.xlsx",      "All_raw",    "DCIS vs Controls (missense)",                          "forest_missense_DCIS_corrected.png",       "Missense (CADD ≥ 20) Variant Associations with DCIS",        "missense",     25,
  "outputs/tables/missense_results_FINAL.xlsx",      "All_raw",    "LCIS vs Controls (missense)",                          "forest_missense_LCIS_corrected.png",       "Missense (CADD ≥ 20) Variant Associations with LCIS",        "missense",     25,
  "outputs/tables/missense_results_FINAL.xlsx",      "All_raw",    "BBD vs Controls (missense)",                           "forest_missense_BBD_corrected.png",        "Missense (CADD ≥ 20) Variant Associations with BBD",         "missense",     25,
  "outputs/tables/missense_results_FINAL.xlsx",      "All_raw",    "Non-proliferative BBD vs Controls (missense)",         "forest_missense_BBD_nonprol_corrected.png","Missense (CADD ≥ 20): Non-proliferative BBD",                "missense",     25,
  "outputs/tables/missense_results_FINAL.xlsx",      "All_raw",    "Proliferative BBD vs Controls (missense)",             "forest_missense_BBD_prol_corrected.png",   "Missense (CADD ≥ 20): Proliferative BBD",                    "missense",     25
)

caption_for <- function(class) if (class == "missense") CAPTION_MISSENSE else CAPTION_TRUNCATING

# Cache each workbook/sheet once.
read_sheet <- local({
  cache <- list()
  function(file, sheet) {
    key <- paste(file, sheet)
    if (is.null(cache[[key]])) cache[[key]] <<- read_excel(file, sheet = sheet)
    cache[[key]]
  }
})

message("Building gene-panel forest figures from saved result tables...")
for (i in seq_len(nrow(SPEC))) {
  s  <- SPEC[i, ]
  df <- read_sheet(s$file, s$sheet) %>% filter(analysis == s$analysis)
  if (nrow(df) == 0) {
    warning(sprintf("No rows for analysis '%s' in %s [%s] — skipped",
                    s$analysis, basename(s$file), s$sheet))
    next
  }
  forest_plot(df, title = s$title,
              file = file.path("outputs/figures", s$out),
              caption = caption_for(s$class), cap = s$cap)
}
message(sprintf("Done: %d forest figures regenerated.", nrow(SPEC)))

# ── Overview heatmap: 9 genes × 6 outcome/variant cells in one figure ────────
message("Building gene × outcome overview heatmap...")
grab <- function(file, sheet, analysis, col_label) {
  read_sheet(file, sheet) %>% filter(analysis == !!analysis) %>%
    transmute(gene, col = col_label, OR, fdr = as.logical(fdr_sig), n_carrier_case)
}
hm_cols <- c("BBD\n(trunc)", "DCIS\n(trunc)", "LCIS\n(trunc)",
             "BBD\n(miss)",  "DCIS\n(miss)",  "LCIS\n(miss)")
hm_df <- bind_rows(
  grab("outputs/tables/BBD_truncating_FINAL.xlsx",         "Raw results", "BBD vs Controls",                 hm_cols[1]),
  grab("outputs/tables/DCIS_LCIS_truncating_results.xlsx", "Raw_results", "DCIS vs Controls (truncating)",   hm_cols[2]),
  grab("outputs/tables/DCIS_LCIS_truncating_results.xlsx", "Raw_results", "LCIS vs Controls (truncating)",   hm_cols[3]),
  grab("outputs/tables/missense_results_FINAL.xlsx",       "All_raw",     "BBD vs Controls (missense)",      hm_cols[4]),
  grab("outputs/tables/missense_results_FINAL.xlsx",       "All_raw",     "DCIS vs Controls (missense)",     hm_cols[5]),
  grab("outputs/tables/missense_results_FINAL.xlsx",       "All_raw",     "LCIS vs Controls (missense)",     hm_cols[6])
)
or_heatmap(hm_df, file = "outputs/figures/heatmap_overview.png",
           title = "Odds ratios across genes × outcomes (truncating and missense)",
           subtitle = paste("* = FDR q < 0.05 (bold).",
                            "† = < 5 carriers among cases: estimate unstable (Firth), colour muted.",
                            "Colour capped at OR 0.25–8."),
           col_levels = hm_cols)

# ── Headline figure: only the FDR-significant DCIS associations ──────────────
message("Building headline summary forest (significant DCIS associations)...")
pick <- function(file, sheet, analysis, genes, suffix) {
  read_sheet(file, sheet) %>%
    filter(analysis == !!analysis, gene %in% genes) %>%
    mutate(gene = paste0(gene, suffix))
}
summary_df <- bind_rows(
  pick("outputs/tables/DCIS_LCIS_truncating_results.xlsx", "Raw_results",
       "DCIS vs Controls (truncating)", c("ATM", "BRCA2", "CHEK2"), " (PTV)"),
  pick("outputs/tables/missense_results_FINAL.xlsx", "All_raw",
       "DCIS vs Controls (missense)", c("CHEK2", "TP53"), " (missense)")
)
summary_order <- c("ATM (PTV)", "BRCA2 (PTV)", "CHEK2 (PTV)",
                   "CHEK2 (missense)", "TP53 (missense)")
forest_plot(summary_df,
            title = "Genes associated with DCIS after FDR correction",
            file  = "outputs/figures/forest_DCIS_significant_summary.png",
            caption = paste(
              "The five DCIS associations significant after BH-FDR (q < 0.05).",
              "PTV = protein-truncating variant; missense = CADD phred ≥ 20. Adjusted for age and study.",
              sep = "\n"),
            genes = summary_order, width = 11)

message("(Hormonal forest is built by objective2_hormonal_analysis.R — different layout.)")
