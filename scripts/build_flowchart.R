# =============================================================================
# Participant flow diagram (STROBE-style) — built from the data, never hand-typed
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# Recomputes every N from the phenotype file and draws the cohort-derivation
# flow with ggplot, so the figure can never drift from the analysis cohorts.
# Output: outputs/figures/participant_flow.png
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2)})
source("scripts/_thesis_theme.R")   # THESIS_PAL, theme helpers, save_fig()
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

ph  <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                  show_col_types = FALSE, na = c("", "NA", "888", "777", "999"))
eur <- ph  %>% filter(ethnicityClass == 1)
ea  <- eur %>% filter(!is.na(ageInt))
cf  <- ea  %>% filter(status == 0)
ins <- ea  %>% filter(status == 2)

N <- list(
  raw      = nrow(ph),
  eur      = nrow(eur),
  excl_eur = nrow(ph)  - nrow(eur),
  elig     = nrow(ea),
  excl_age = nrow(eur) - nrow(ea),
  cf       = nrow(cf),
  insitu   = nrow(ins),
  inv      = sum(ea$status == 1, na.rm = TRUE),
  dcis     = sum(ins$MorphologygroupIndex_corr == "Ductal",  na.rm = TRUE),
  lcis     = sum(ins$MorphologygroupIndex_corr == "Lobular", na.rm = TRUE),
  bbd_case = sum(cf$BBD_history == 1, na.rm = TRUE),
  bbd_ctrl = sum(cf$BBD_history == 0, na.rm = TRUE)
)
N$insitu_other <- N$insitu - N$dcis - N$lcis
N$bbd_unknown  <- N$cf - N$bbd_case - N$bbd_ctrl
f <- function(x) formatC(x, big.mark = ",", format = "d")

NAVY <- unname(THESIS_PAL["baseline"]); GREY <- "grey45"
# Boxes: id, x, y (centre), width, height, label, kind {step, exclude, group, outcome}
B <- tribble(
  ~x,  ~y,  ~w,  ~h, ~label,                                                            ~kind,
  31,  93,  34,  9, sprintf("All BRIDGES participants\n(phenotype file v17)\nn = %s", f(N$raw)),                "step",
  31,  76,  34,  9, sprintf("European ancestry\n(ethnicityClass = 1)\nn = %s", f(N$eur)),                       "step",
  31,  59,  34,  9, sprintf("Analysis-eligible\n(non-missing age at interview)\nn = %s", f(N$elig)),            "step",
  80,  84.5,30,  7, sprintf("Excluded: non-European\nn = %s", f(N$excl_eur)),                                   "exclude",
  80,  67.5,30,  7, sprintf("Excluded: missing age\nn = %s", f(N$excl_age)),                                    "exclude",
  14,  41,  24,  9, sprintf("Cancer-free\n(status = 0)\nn = %s", f(N$cf)),                                       "group",
  46,  41,  24,  9, sprintf("In-situ carcinoma\n(status = 2)\nn = %s", f(N$insitu)),                            "group",
  80,  41,  28, 10, sprintf("Invasive BC\n(status = 1)\nn = %s\n[positive control]", f(N$inv)),                 "group",
  8,   18,  18,  9, sprintf("BBD cases\n(BBD history +)\nn = %s", f(N$bbd_case)),                               "outcome",
  28,  18,  20,  9, sprintf("BBD controls\n(BBD history −)\nn = %s", f(N$bbd_ctrl)),                            "outcome",
  41,  18,  17,  9, sprintf("DCIS\n(Ductal)\nn = %s", f(N$dcis)),                                                "outcome",
  60,  18,  17,  9, sprintf("LCIS\n(Lobular)\nn = %s", f(N$lcis)),                                               "outcome"
)
fill_for   <- c(step = "white", exclude = "grey94", group = "#D6E0F0", outcome = NAVY)
border_for <- c(step = NAVY,    exclude = GREY,     group = NAVY,      outcome = NAVY)
text_for   <- c(step = "black", exclude = GREY,     group = NAVY,      outcome = "white")
B <- B %>% mutate(fill = fill_for[kind], border = border_for[kind],
                  tcol = text_for[kind], lty = ifelse(kind == "exclude", "22", "solid"))

# Arrows: x,y -> xend,yend
seg <- function(x, y, xe, ye, dashed = FALSE) tibble(x, y, xe, ye, dashed)
A <- bind_rows(
  seg(31, 93 - 4.5, 31, 76 + 4.5),                  # raw -> european
  seg(31, 76 - 4.5, 31, 59 + 4.5),                  # european -> eligible
  seg(31, 84.5, 65, 84.5, TRUE),                    # -> exclude non-eur
  seg(31, 67.5, 65, 67.5, TRUE),                    # -> exclude age
  seg(31, 59 - 4.5, 14, 41 + 4.5),                  # eligible -> cancer-free
  seg(31, 59 - 4.5, 46, 41 + 4.5),                  # eligible -> in-situ
  seg(31, 59 - 4.5, 80, 41 + 5),                    # eligible -> invasive
  seg(14, 41 - 4.5, 8,  18 + 4.5),                  # cancer-free -> bbd cases
  seg(14, 41 - 4.5, 28, 18 + 4.5),                  # cancer-free -> bbd controls
  seg(46, 41 - 4.5, 41, 18 + 4.5),                  # in-situ -> dcis
  seg(46, 41 - 4.5, 60, 18 + 4.5)                   # in-situ -> lcis
)

p <- ggplot() +
  geom_segment(data = A, aes(x = x, y = y, xend = xe, yend = ye,
                             linetype = dashed),
               arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
               colour = GREY, linewidth = 0.5, show.legend = FALSE) +
  geom_rect(data = B, aes(xmin = x - w/2, xmax = x + w/2,
                          ymin = y - h/2, ymax = y + h/2),
            fill = B$fill, colour = B$border, linewidth = 0.6) +
  geom_text(data = B, aes(x = x, y = y, label = label),
            colour = B$tcol, size = 2.85, lineheight = 0.95, fontface = "plain") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22")) +
  annotate("text", x = 14, y = 41 - 6.2, hjust = 0.5, size = 2.4, colour = GREY,
           label = sprintf("(also the control pool for DCIS & LCIS; BBD history unknown excluded, n = %s)",
                           f(N$bbd_unknown))) +
  annotate("text", x = 50.5, y = 18 - 6.2, hjust = 0.5, size = 2.4, colour = GREY,
           label = sprintf("(other / mixed in-situ morphology excluded, n = %s)", f(N$insitu_other))) +
  labs(title = "Derivation of the analysis cohorts",
       subtitle = "BRIDGES phenotype file v17; European-ancestry, age-complete participants.") +
  coord_cartesian(xlim = c(-2, 96), ylim = c(8, 99)) +
  theme_void(base_size = 12, base_family = THESIS_FONT) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 8.5, colour = GREY, hjust = 0),
        plot.margin   = margin(10, 10, 10, 10))

save_fig(p, "outputs/figures/participant_flow.png", width = 11, height = 8)
