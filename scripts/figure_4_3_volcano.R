# =============================================================================
# Figure 4.3 — BFDP-encoded volcano plot of DCIS gene associations
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil Population Health Sciences — 2025/26
#
# PURPOSE
#   A single figure that places every DCIS gene-level association (nine genes
#   x truncating/missense = 18 points) in effect-size / significance space and
#   simultaneously encodes its Bayesian credibility. Four quantities are shown:
#     x  = log2(odds ratio)
#     y  = -log10(p-value)
#     colour = BFDP (continuous; darker = less credible), threshold at 0.05
#     shape  = variant class (circle = truncating, triangle = missense)
#   Point size is fixed; carrier count is available in Table 4.2/4.3/4.3a.
#
#   Reads BFDP from outputs/tables/BFDP_analysis_DCIS.csv (produced by
#   bfdp_analysis.R) so the figure and the table can never disagree. Appearance
#   follows scripts/_thesis_theme.R; this script never alters a reported value.
#
# OUTPUT: outputs/figures/figure_4_3_volcano.png (+ vector .pdf)
# =============================================================================

suppressMessages({ library(ggplot2); library(dplyr) })
source("scripts/_thesis_theme.R")

# -- 1. DATA ------------------------------------------------------------------
d <- read.csv("outputs/tables/BFDP_analysis_DCIS.csv", check.names = FALSE) %>%
  rename(p = p, BFDP = BFDP) %>%
  mutate(
    log2OR   = log2(OR),
    neglog10p = -log10(p),
    variant  = factor(Variant, levels = c("truncating", "missense"),
                       labels = c("Truncating", "Missense")),
    sig      = `BH/FDR` == "BH-sig",
    # BFDP is reported rounded to 3 dp; ATM rounds to 0.000. Floor at 1e-4 so
    # the log-scale colour maps it to the brightest (most-credible) end rather
    # than to NA/grey. This affects the COLOUR channel only, not any value.
    BFDP_col = pmax(BFDP, 1e-4)
  )

# -- 2. PLOT ------------------------------------------------------------------
p <- ggplot(d, aes(log2OR, neglog10p)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = THESIS_PAL["rule"],
             linewidth = 0.4) +
  geom_point(aes(colour = BFDP_col, shape = variant), size = 3.2, stroke = 0.4) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(sig | BFDP > 0.4,
                       paste0(Gene, " (", ifelse(is.na(`Carrier cases`), "", `Carrier cases`),
                              ifelse(is.na(`Carrier cases`), "", " carr."), ")"), Gene)),
    size = 2.9, family = THESIS_FONT, colour = "grey20",
    max.overlaps = 20, box.padding = 0.4, min.segment.length = 0.2, seed = 1) +
  scale_colour_viridis_c(
    option = "viridis", direction = -1, trans = "log10",
    breaks = c(0.001, 0.01, 0.05, 0.2, 0.5),
    labels = c("0.001", "0.01", "0.05", "0.2", "0.5"),
    name = "BFDP") +
  scale_shape_manual(values = c("Truncating" = 16, "Missense" = 17), name = NULL) +
  labs(
    x = expression(log[2]~"(odds ratio)"),
    y = expression(-log[10]~"(" * italic(p) * "-value)"),
    title = "DCIS gene associations: effect size, significance and credibility",
    caption = paste0("Colour = Wakefield BFDP (darker = less credible; ",
                     "threshold 0.05). Shape = variant class. TP53 truncating ",
                     "lies far right (OR 5.88) but is\ndark (BFDP 0.44) — one carrier, ",
                     "not credible. Five associations are BH-significant (Table 4.3a).")) +
  theme_thesis(grid = "none") +
  theme(legend.position = "right", legend.title = element_text())

# -- 3. SAVE ------------------------------------------------------------------
if (!requireNamespace("ggrepel", quietly = TRUE))
  stop("ggrepel required; install with manage_packages r-ggrepel")
save_fig(p, "outputs/figures/figure_4_3_volcano.png", width = 7.4, height = 5.2)
cat("Figure 4.3 written. 18 points; five BH-significant:",
    paste(d$Gene[d$sig], collapse = ", "), "\n")
