# =============================================================================
# Stage-of-progression slopegraph: where in the benign -> in-situ sequence does
# each truncating-variant effect emerge?
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# The headline finding is not just "ATM/CHEK2/BRCA2 are associated with DCIS" but
# WHERE on the histological sequence the signal appears. This figure places the
# odds ratio for each gene at two adequately powered stages — benign breast
# disease (cancer-free women with a BBD history) and DCIS — on one log axis, so a
# reader sees the pattern directly: ATM and CHEK2 are null at the benign stage
# and switch on at DCIS, whereas BRCA2 is already elevated at the benign stage.
#
# Only stages with a usable number of carriers are plotted. The BBD histological
# SUBTYPES (non-proliferative / proliferative) and every LCIS cell rest on 0-3
# carriers, so their Firth point estimates are unstable and are deliberately not
# drawn here (they remain in the full result tables). No model is re-fitted; the
# numbers are read from the saved, verified result workbooks.
#
# Run:  Rscript scripts/build_gradient.R
# =============================================================================

suppressPackageStartupMessages({library(readxl); library(dplyr); library(ggrepel)})
source("scripts/_thesis_theme.R")
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

GENES <- c("ATM", "CHEK2", "BRCA2")           # the FDR-significant PTV genes
STAGE_LEVELS <- c("Benign breast disease", "DCIS")

# One row per (gene, stage): read straight from the verified result tables.
bbd  <- read_excel("outputs/tables/BBD_truncating_FINAL.xlsx",         sheet = "Raw results")
dcis <- read_excel("outputs/tables/DCIS_LCIS_truncating_results.xlsx", sheet = "Raw_results")

pull <- function(tbl, analysis, stage) {
  tbl %>%
    filter(analysis == !!analysis, gene %in% GENES) %>%
    transmute(gene, stage = stage, OR, CI_low, CI_high,
              fdr = as.logical(fdr_sig), n_carrier_case)
}

grad <- bind_rows(
  pull(bbd,  "BBD vs Controls",               STAGE_LEVELS[1]),
  pull(dcis, "DCIS vs Controls (truncating)", STAGE_LEVELS[2])
) %>%
  mutate(
    gene   = factor(gene, levels = GENES),
    stage  = factor(stage, levels = STAGE_LEVELS),
    stage_x = as.integer(stage),
    sig_lab = if_else(fdr, "FDR q < 0.05", "Not significant")
  )

GENE_COLOURS <- setNames(THESIS_GROUPS, GENES)   # orange / dark-blue / light-blue
SIG_SHAPES   <- c("FDR q < 0.05" = 18, "Not significant" = 5)  # filled vs hollow diamond

# Right-hand direct labels (slopegraph convention: no colour legend needed).
end_labels <- grad %>% filter(stage == "DCIS")

p <- ggplot(grad, aes(x = stage_x, y = OR, colour = gene, group = gene)) +
  geom_hline(yintercept = 1, linetype = "longdash",
             colour = unname(THESIS_PAL["rule"]), linewidth = 0.5) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.05, linewidth = 0.7) +
  geom_point(aes(shape = sig_lab, fill = gene), size = 4.2, stroke = 1) +
  geom_text_repel(data = end_labels,
                  aes(label = sprintf("%s  %.2f", gene, OR)),
                  hjust = 0, direction = "y", nudge_x = 0.10,
                  segment.colour = "grey75", segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.35,
                  size = 3.6, fontface = "bold") +
  scale_x_continuous(breaks = c(1, 2), labels = STAGE_LEVELS,
                     limits = c(0.85, 2.75)) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  scale_colour_manual(values = GENE_COLOURS, guide = "none") +
  scale_fill_manual(values = GENE_COLOURS, guide = "none") +
  scale_shape_manual(values = SIG_SHAPES, name = NULL) +
  labs(
    title = "Where the truncating-variant effect emerges along the progression",
    subtitle = "Odds ratio for each gene at two adequately powered stages (log scale). Filled = FDR q < 0.05.",
    x = NULL, y = "Odds ratio vs cancer-free controls",
    caption = paste(
      "Firth's penalised logistic regression, adjusted for age and study; 95% profile-likelihood CIs.",
      "ATM and CHEK2 are null at the benign stage and only reach significance at DCIS; BRCA2 is elevated at both.",
      "BBD histological subtypes and LCIS rest on 0-3 carriers (unstable) and are not plotted.",
      sep = "\n")
  ) +
  theme_thesis(base_size = 13, grid = "y") +
  theme(axis.text.x = element_text(face = "bold", size = 12))

save_fig(p, "outputs/figures/gradient_stage_truncating.png", width = 9, height = 6.5)
message("Done: stage-of-progression slopegraph written.")
