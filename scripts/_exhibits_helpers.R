# =============================================================================
# Shared exhibit helpers — one forest-plot function for every gene-panel figure
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# WHY THIS FILE EXISTS
#   The DCIS, BBD and missense analysis scripts each carried a near-identical
#   ~100-line forest-plot function that had drifted apart in style (theme_bw vs
#   theme_classic, different CI capping, different text-column positions). This
#   one function gives every panel figure a single consistent appearance, so a
#   styling change is made in one place rather than three.
#
#   Crucially, figures are drawn from the SAVED, verified result tables
#   (outputs/tables/*.xlsx) rather than by re-fitting any model, so
#   regenerating the exhibits cannot change a single reported number.
#
# CONTRACT
#   forest_plot(df, title, file, caption, cap, width, genes)
#     df : data frame with columns
#          gene, OR, CI_low, CI_high, fdr_sig, n_carrier_case, n_carrier_ctrl
#   Side effect: writes a 300-dpi PNG to `file`. Returns the ggplot invisibly.
# =============================================================================

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("scripts/_thesis_theme.R")   # THESIS_PAL, theme_thesis(), save_fig()

# Panel gene order (top-to-bottom on the plot is the reverse of this).
GENES_PANEL <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
                 "BARD1", "RAD51C", "RAD51D", "TP53")

CAPTION_TRUNCATING <- paste(
  "Firth's penalised logistic regression (profile penalised-likelihood CIs). Adjusted for age and study.",
  "BH-FDR across 9 genes. Red = FDR q < 0.05; arrow = CI truncated at panel edge.",
  sep = "\n")

CAPTION_MISSENSE <- paste(
  "Firth's penalised logistic regression (profile penalised-likelihood CIs). Carriers: CADD phred ≥ 20 (top 1% deleterious).",
  "Adjusted for age and study. BH-FDR across 9 genes. Red = FDR q < 0.05; arrow = CI truncated.",
  sep = "\n")

# Colour/shape conventions shared by every panel figure (from the house palette).
.FOREST_COLOURS <- c("FDR significant" = unname(THESIS_PAL["signal"]),
                     "Not significant" = unname(THESIS_PAL["baseline"]),
                     "Not estimable"   = unname(THESIS_PAL["muted"]))
.FOREST_SHAPES  <- c("FDR significant" = 18,
                     "Not significant" = 18,
                     "Not estimable"   = 1)

forest_plot <- function(df, title, file,
                        caption = CAPTION_TRUNCATING,
                        cap     = 25,
                        width   = 12,
                        genes   = GENES_PANEL) {

  sig <- !is.na(df$fdr_sig) & as.logical(df$fdr_sig)

  plot_df <- df %>%
    mutate(
      gene_f       = factor(gene, levels = rev(genes)),
      colour_group = case_when(is.na(OR) ~ "Not estimable",
                               sig       ~ "FDR significant",
                               TRUE      ~ "Not significant"),
      or_label      = if_else(is.na(OR), "Not estimable",
                              sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)),
      carrier_label = sprintf("%d / %d", n_carrier_case, n_carrier_ctrl)
    )

  # Adaptive panel width: cap the x-axis so a few artifact CIs (e.g. LCIS
  # TP53/BARD1, upper bound > 200) do not compress the log axis. The full CI is
  # still printed in the text column; a small arrow flags any bar clipped here.
  x_min <- min(c(plot_df$CI_low, 0.3), na.rm = TRUE) * 0.75
  x_max <- min(max(c(plot_df$CI_high, 5.0), na.rm = TRUE), cap) * 1.20
  plot_df <- plot_df %>%
    mutate(CI_high_draw = pmin(CI_high, x_max),
           clipped      = !is.na(CI_high) & CI_high > x_max)

  x_carrier <- x_max * 1.30
  x_or      <- x_max * 2.40

  p <- ggplot(plot_df, aes(y = gene_f)) +
    geom_vline(xintercept = 1, linetype = "longdash",
               colour = unname(THESIS_PAL["rule"]), linewidth = 0.5) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high_draw, colour = colour_group),
                  width = 0.28, linewidth = 0.8, na.rm = TRUE, orientation = "y") +
    geom_text(data = ~ filter(.x, clipped),
              aes(x = x_max, label = "→"),
              hjust = 0, vjust = 0.35, size = 3.2,
              colour = unname(THESIS_PAL["baseline"]), na.rm = TRUE) +
    geom_point(aes(x = OR, colour = colour_group, shape = colour_group),
               size = 4, na.rm = TRUE) +
    geom_text(aes(x = x_carrier, label = carrier_label),
              hjust = 0.5, size = 3, colour = "grey25") +
    geom_text(aes(x = x_or, label = or_label),
              hjust = 0, size = 3, colour = "grey15") +
    annotate("text", x = x_carrier, y = length(genes) + 0.7,
             label = "Cases / Ctrl", size = 3, fontface = "bold", hjust = 0.5) +
    annotate("text", x = x_or, y = length(genes) + 0.7,
             label = "OR (95% CI)", size = 3, fontface = "bold", hjust = 0) +
    scale_x_log10(breaks = c(0.3, 0.5, 1, 2, 5, 10, 20),
                  labels = c("0.3", "0.5", "1", "2", "5", "10", "20"),
                  limits = c(x_min, x_or * 1.9)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
    scale_colour_manual(values = .FOREST_COLOURS, name = NULL) +
    scale_shape_manual(values = .FOREST_SHAPES, name = NULL) +
    labs(title = title, x = "Odds Ratio (log scale, 95% profile-likelihood CI)", y = NULL,
         caption = caption) +
    theme_thesis(base_size = 13, grid = "y") +
    theme(axis.text.y = element_text(face = "bold", size = 11))

  save_fig(p, file, width = width, height = 6.5)
  invisible(p)
}

# -----------------------------------------------------------------------------
# or_heatmap(): compact landscape of every gene × outcome odds ratio in one
# figure. `df` is long with columns gene, col (outcome label), OR, fdr (logical),
# and optionally n_carrier_case.
# Colour is log2(OR), capped to OR 0.25–8 so artifact CIs don't blow out the
# scale; FDR-significant cells are bold and marked with an asterisk; cells with
# no carriers among cases are labelled "NE" (not estimable).
# When n_carrier_case is supplied, cells with fewer than `min_carrier_case`
# carriers among cases carry an unreliable Firth point estimate (e.g. LCIS
# BARD1 OR 19.81 from 0 case carriers). These are marked with a dagger and their
# fill is muted, so a near-empty cell no longer reads as a strong colour signal.
# -----------------------------------------------------------------------------
or_heatmap <- function(df, file,
                       title    = "Odds ratios across genes and outcomes",
                       subtitle = "* = FDR q < 0.05 (bold). † = < 5 carriers among cases (estimate unstable, muted). NE = not estimable. Colour capped at OR 0.25–8.",
                       genes      = GENES_PANEL,
                       col_levels = NULL,
                       width = 9, height = 6,
                       min_carrier_case = 5) {
  if (is.null(col_levels)) col_levels <- unique(df$col)
  has_carrier <- "n_carrier_case" %in% names(df)
  hm <- df %>%
    mutate(gene = factor(gene, levels = rev(genes)),
           col  = factor(col, levels = col_levels),
           bold = !is.na(fdr) & as.logical(fdr),
           unstable = if (has_carrier) !is.na(n_carrier_case) & n_carrier_case < min_carrier_case else FALSE,
           l2   = pmax(pmin(log2(OR), 3), -2),
           lab  = ifelse(is.na(OR), "NE",
                         sprintf("%.2f%s%s", OR, ifelse(bold, "*", ""), ifelse(unstable, "†", ""))))

  p <- ggplot(hm, aes(col, gene, fill = l2)) +
    geom_tile(aes(alpha = unstable), colour = "white", linewidth = 1.1) +
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.25), guide = "none") +
    geom_text(aes(label = lab, fontface = ifelse(bold, "bold", "plain")),
              colour = "grey10", size = 3.4) +
    scale_fill_gradient2(low = unname(THESIS_DIVERGING["low"]),
                         mid = unname(THESIS_DIVERGING["mid"]),
                         high = unname(THESIS_DIVERGING["high"]),
                         midpoint = 0, na.value = "grey80", name = "Odds ratio",
                         breaks = c(-2, 0, 3), labels = c("0.25", "1", "8")) +
    scale_x_discrete(position = "top") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_thesis(base_size = 13, grid = "none") +
    theme(axis.text.y     = element_text(face = "bold"),
          axis.ticks      = element_blank(),
          legend.position = "right",
          legend.title    = element_text(size = 11, colour = "grey15"))

  save_fig(p, file, width = width, height = height)
  invisible(p)
}
