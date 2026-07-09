# =============================================================================
# Control-Group Sensitivity — DCIS genetic findings with BBD-free controls
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE (addresses supervisor point #8: "choice of control group")
#   There is a cross-objective inconsistency:
#     - BBD analysis controls  = cancer-free women WITHOUT BBD history
#                                (status==0 & BBD_history==0)
#     - DCIS/LCIS controls      = FULL cancer-free pool (status==0), which
#                                INCLUDES 1,765 confirmed-BBD women.
#   If BBD shares genetic susceptibility with DCIS, including BBD women in the
#   DCIS control group would bias DCIS odds ratios toward the null (controls
#   "contaminated" with susceptibility carriers). This script re-estimates the
#   five FDR-significant DCIS findings under two control definitions:
#     (P) Primary controls  : status==0                       [as reported]
#     (X) BBD-free controls  : status==0 & NOT confirmed BBD   [removes 1,765]
#   Women with UNKNOWN BBD history (BBD_history == NA, n=19,944) are retained
#   in (X): the test isolates the effect of the 1,765 KNOWN-BBD controls, and
#   does not confound with the differential collection of BBD history by study.
#
#   Stable ORs => the control-pool inconsistency is immaterial. A systematic
#   strengthening under (X) => the primary DCIS ORs are mildly diluted by
#   BBD contamination of controls (conservative bias, reassuring direction).
#
# Built on sensitivity_familial_population.R machinery (same Firth estimator).
# DATA GOVERNANCE: BCAC/BRIDGES approved secondary use; non-disclosive output.
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(logistf)
  library(writexl)
})

PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
MISS_FILE  <- "concept_807_zhang_bridges_missense.csv"
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

GENES <- c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM",
           "BARD1", "RAD51C", "RAD51D", "TP53")

# Five FDR-significant DCIS findings to defend
TARGETS <- tribble(
  ~vtype,        ~gene,
  "truncating",  "ATM",
  "truncating",  "CHEK2",
  "truncating",  "BRCA2",
  "missense",    "CHEK2",
  "missense",    "TP53"
)

pheno <- suppressWarnings(read_delim(PHENO_FILE, delim = "\t",
  show_col_types = FALSE, na = c("", "NA", "888", "777", "999")))
trunc <- read_csv(TRUNC_FILE, show_col_types = FALSE, na = c("", "NA")) %>% select(-any_of("...1"))
miss  <- read_csv(MISS_FILE,  show_col_types = FALSE, na = c("", "NA")) %>% select(-any_of("...1"))

build_carriers <- function(raw, pattern, suffix) {
  cols <- map_chr(GENES, function(g) {
    hit <- names(raw)[grepl(paste0("^", g, pattern), names(raw), ignore.case = TRUE)]
    if (length(hit) > 0) hit[1] else NA_character_
  })
  keep <- !is.na(cols)
  raw %>% select(BRIDGES_ID, all_of(cols[keep])) %>%
    rename(!!!setNames(cols[keep], paste0(GENES[keep], suffix))) %>%
    mutate(across(-BRIDGES_ID, ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))
}
carriers_t <- build_carriers(trunc, "_truncating$",        "_t")
carriers_m <- build_carriers(miss,  "_CADD\\.phred\\.01$", "_m")

dat <- pheno %>%
  filter(ethnicityClass == 1, !is.na(ageInt)) %>%
  left_join(carriers_t, by = "BRIDGES_ID") %>%
  left_join(carriers_m, by = "BRIDGES_ID") %>%
  mutate(across(c(ends_with("_t"), ends_with("_m")), ~ replace_na(.x, 0L)))

dcis            <- dat %>% filter(status == 2, MorphologygroupIndex_corr == "Ductal")
controls_full   <- dat %>% filter(status == 0)                                   # (P)
controls_bbdfree<- dat %>% filter(status == 0, is.na(BBD_history) | BBD_history == 0)  # (X)

message(sprintf("Controls (P) full: %d | (X) BBD-free: %d | removed confirmed-BBD: %d",
                nrow(controls_full), nrow(controls_bbdfree),
                nrow(controls_full) - nrow(controls_bbdfree)))

fit_ctrl <- function(gene, vtype, ctrl) {
  carrier_col <- paste0(gene, if (vtype == "missense") "_m" else "_t")
  d <- bind_rows(ctrl %>% mutate(outcome = 0L),
                 dcis %>% mutate(outcome = 1L)) %>%
    mutate(carrier = .data[[carrier_col]], study = droplevels(factor(study)))
  fit <- logistf(outcome ~ carrier + ageInt + study, data = d, firth = TRUE, pl = FALSE)
  i <- which(names(coef(fit)) == "carrier")
  tibble(OR = exp(coef(fit)[i]), CI_low = exp(fit$ci.lower[i]),
         CI_high = exp(fit$ci.upper[i]), p = fit$prob[i],
         carr_ctrl = sum(d$carrier == 1 & d$outcome == 0))
}

res <- TARGETS %>%
  mutate(P = pmap(list(gene, vtype), ~ fit_ctrl(..1, ..2, controls_full)),
         X = pmap(list(gene, vtype), ~ fit_ctrl(..1, ..2, controls_bbdfree))) %>%
  mutate(
    `OR (P) primary`   = map_chr(P, ~ sprintf("%.2f (%.2f–%.2f)", .$OR, .$CI_low, .$CI_high)),
    `p (P)`            = map_chr(P, ~ sprintf("%.3f", .$p)),
    `OR (X) BBD-free`  = map_chr(X, ~ sprintf("%.2f (%.2f–%.2f)", .$OR, .$CI_low, .$CI_high)),
    `p (X)`            = map_chr(X, ~ sprintf("%.3f", .$p)),
    rel_shift = map2_dbl(P, X, ~ (.y$OR - .x$OR) / .x$OR),
    direction = if_else(rel_shift >= 0, "stronger (as expected if dilution)", "weaker"),
    `OR shift` = sprintf("%+.1f%%", 100 * rel_shift)
  ) %>%
  select(vtype, gene, `OR (P) primary`, `p (P)`,
         `OR (X) BBD-free`, `p (X)`, `OR shift`, direction)

cat("\n", strrep("=", 96), "\n",
    "CONTROL-GROUP SENSITIVITY — 5 FDR-significant DCIS findings, full vs BBD-free controls\n",
    strrep("=", 96), "\n", sep = "")
as.data.frame(res) %>% print(right = FALSE)

write_xlsx(list("Control_sensitivity" = as.data.frame(res)),
           "outputs/tables/control_group_sensitivity.xlsx")
message("\nExported: outputs/tables/control_group_sensitivity.xlsx")
