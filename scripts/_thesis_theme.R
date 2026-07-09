# =============================================================================
# Shared house style for every dissertation figure
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# WHY THIS FILE EXISTS
#   The figure set was drawn by several scripts that each picked their own
#   ggplot theme (theme_bw, theme_grey, theme_minimal) and their own colours.
#   A single house theme + one colour dictionary gives the whole thesis one
#   visual voice, which is exactly what a reader scores on first glance.
#
#   This file changes ONLY appearance. It never touches a reported number:
#   figures are still drawn from the saved, verified result tables.
#
# WHAT IT PROVIDES
#   THESIS_FONT  : the sans-serif family used everywhere (Arial)
#   THESIS_PAL   : a named colour dictionary with one fixed meaning per colour
#   theme_thesis(): the Nature-leaning base theme (white, faint guides, Arial)
#   save_fig()   : write a high-DPI PNG preview AND an editable vector PDF
# =============================================================================

suppressPackageStartupMessages({library(ggplot2)})

# Sans-serif family. Arial is present on macOS and embeds cleanly in cairo_pdf;
# fall back to the device default if a machine lacks it.
THESIS_FONT <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else ""

# Colour dictionary — ONE meaning per colour, applied across every figure.
#   signal   red   : a highlighted / statistically significant result
#   baseline blue  : the non-significant / reference category
#   muted    grey  : not estimable or deliberately de-emphasised
#   accent   orange + accent_alt light-blue : extra subgroup encodings
#                    (e.g. the three study-inclusion specs in sensitivity plots)
#   rule     grey  : the OR = 1 reference line
# NOTE: a diverging heatmap legitimately uses red=high / blue=low as a
# CONTINUOUS effect-size scale; that is a different channel from the categorical
# "red = significant" above and the two never appear in the same panel.
THESIS_PAL <- c(
  signal     = "#B03A2E",
  baseline   = "#1A5276",
  muted      = "grey70",
  accent     = "#CA6B27",
  accent_alt = "#5499C7",
  rule       = "grey55"
)

# Diverging fill endpoints for heatmaps (low OR -> high OR through neutral).
THESIS_DIVERGING <- c(low = "#1A5276", mid = "grey95", high = "#B03A2E")

# Qualitative colours for 2-3 NON-significance subgroups (e.g. study-inclusion
# specifications). Deliberately excludes the significance-red so the meaning of
# red stays "significant / highlighted" everywhere else.
THESIS_GROUPS <- unname(THESIS_PAL[c("baseline", "accent", "accent_alt")])

# Make every in-plot text/label use the house font too (geoms ignore the theme's
# base_family by default). Sourcing this file is enough.
update_geom_defaults("text",  list(family = THESIS_FONT))
update_geom_defaults("label", list(family = THESIS_FONT))

# Nature-leaning base theme: white panel, no box, faint guide lines only on the
# requested axis, outward ticks, Arial, left-aligned bold title.
#   grid = "y"   horizontal guides only (forest plots, dot plots)
#   grid = "x"   vertical guides only
#   grid = "none" no guides (heatmaps, schematics)
theme_thesis <- function(base_size = 12, base_family = THESIS_FONT, grid = "y") {
  th <- theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(colour = "grey15"),
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 2, colour = "grey35", hjust = 0,
                                      margin = margin(b = 6)),
      plot.caption     = element_text(size = base_size - 3.5, colour = "grey45", hjust = 0,
                                      margin = margin(t = 8)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title       = element_text(colour = "grey15", size = base_size),
      axis.text        = element_text(colour = "grey25"),
      axis.ticks       = element_line(colour = "grey70", linewidth = 0.3),
      axis.ticks.length = unit(2.5, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.key.size  = unit(11, "pt"),
      legend.text      = element_text(size = base_size - 1),
      strip.text       = element_text(face = "bold", size = base_size - 1, colour = "grey15"),
      strip.background = element_blank()
    )
  if (grid == "y")    th <- th + theme(panel.grid.major.x = element_blank())
  if (grid == "x")    th <- th + theme(panel.grid.major.y = element_blank())
  if (grid == "none") th <- th + theme(panel.grid.major   = element_blank())
  th
}

# Write a figure as BOTH a high-DPI PNG (for Word / quick view) and an editable
# vector PDF (for the submission PDF — crisp text and lines at any size).
# `file` is the .png path; the .pdf is written alongside with the same stem.
save_fig <- function(plot, file, width, height, dpi = 320) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(file, width = width, height = height, units = "in",
                res = dpi, background = "white")
  print(plot); dev.off()
  pdf_file <- sub("\\.png$", ".pdf", file)
  grDevices::cairo_pdf(pdf_file, width = width, height = height, family = THESIS_FONT)
  print(plot); dev.off()
  message(sprintf("  %s  (+ vector %s)", basename(file), basename(pdf_file)))
  invisible(plot)
}
