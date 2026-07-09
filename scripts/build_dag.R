# =============================================================================
# Collider-bias DAG for the benign breast disease (BBD) genetic analysis
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# Visualises why the BBD germline-association estimates are biased. The BBD
# analysis is restricted to cancer-free women (status == 0). Panel-gene carriers
# and BBD history are both causes of breast cancer, so breast cancer is a
# COLLIDER on the path  gene -> breast cancer <- BBD. Conditioning the sample on
# cancer-free status (a function of breast cancer) opens a spurious, non-causal
# association between carrier status and BBD — carriers who progressed to cancer
# are selectively removed — which produces the observed below-unity odds ratios.
#
# Output: outputs/figures/collider_dag.png  (a schematic; no data are read)
# =============================================================================

suppressPackageStartupMessages({library(ggplot2)})
source("scripts/_thesis_theme.R")   # THESIS_PAL, theme helpers, save_fig()
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

NAVY <- unname(THESIS_PAL["baseline"]); RED <- unname(THESIS_PAL["signal"]); GREY <- "grey45"

nodes <- data.frame(
  x     = c(1.0, 5.0, 3.0),
  y     = c(3.0, 3.0, 1.0),
  label = c("Panel gene\n(carrier)", "BBD history", "Breast cancer"),
  id    = c("G", "B", "C"))

# Causal edges (solid, black) -> the collider.
causal <- data.frame(x = c(1.25, 4.75), y = c(2.75, 2.75),
                      xe = c(2.7, 3.3),  ye = c(1.25, 1.25))

p <- ggplot() +
  # spurious, non-causal association induced by conditioning on the collider
  geom_curve(aes(x = 1.45, y = 3.2, xend = 4.55, yend = 3.2),
             curvature = -0.35, linetype = "22", colour = RED, linewidth = 0.7) +
  annotate("text", x = 3.0, y = 4.3, colour = RED, size = 3.2, fontface = "bold",
           label = "spurious association induced by conditioning") +
  # causal arrows into the collider
  geom_segment(data = causal, aes(x = x, y = y, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "black", linewidth = 0.7) +
  # selection box around the collider (sample conditioned on cancer-free)
  annotate("rect", xmin = 2.05, xmax = 3.95, ymin = 0.6, ymax = 1.4,
           fill = NA, colour = RED, linewidth = 0.8, linetype = "solid") +
  annotate("text", x = 3.0, y = 0.35, colour = RED, size = 2.9,
           label = "[ conditioned: sample restricted to cancer-free women ]") +
  # nodes
  geom_label(data = nodes, aes(x = x, y = y, label = label),
             fill = "white", colour = NAVY, fontface = "bold",
             label.size = 0.6, size = 3.4, label.padding = unit(0.4, "lines")) +
  annotate("text", x = 3.0, y = -0.35, colour = GREY, size = 2.9, hjust = 0.5,
           label = paste("Breast cancer is a collider on  gene → breast cancer ← BBD.  Restricting to cancer-free women conditions on a function of",
                         "the collider, opening the dashed path: carriers who progressed to cancer leave the cancer-free group, depleting",
                         "carriers among BBD-history cases and biasing the gene–BBD odds ratio below 1.", sep = "\n")) +
  labs(title = "Collider bias in the cancer-free BBD analysis") +
  coord_cartesian(xlim = c(0, 6), ylim = c(-1.1, 4.7)) +
  theme_void(base_size = 12, base_family = THESIS_FONT) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.margin = margin(10, 14, 10, 14))

save_fig(p, "outputs/figures/collider_dag.png", width = 9, height = 6.2)
