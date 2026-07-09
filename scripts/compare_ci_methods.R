# =============================================================================
# CI-method check — Wald vs profile-likelihood for the FDR-significant DCIS genes
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — University of Cambridge MPhil
#
# PURPOSE
#   The DCIS/LCIS truncating analysis uses Wald CIs (pl = FALSE) because
#   profile-likelihood (pl = TRUE) profiling hangs on the LCIS genes that have
#   0 carriers among cases. This script confirms that, for the three genes that
#   actually drive the DCIS findings (ATM, CHEK2, BRCA2 — all with ample
#   carriers), Wald and profile-likelihood CIs are materially identical, so the
#   Wald choice does not affect any reported conclusion.
#
# OUTPUT
#   outputs/tables/ci_method_check.xlsx  (Wald vs PL side by side)
# =============================================================================

required <- c("tidyverse", "logistf", "writexl")
for (pkg in required)
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
library(tidyverse); library(logistf); library(writexl)

PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
SIG_GENES  <- c("ATM", "CHEK2", "BRCA2")   # the FDR-significant DCIS truncating genes

pheno <- read_delim(PHENO_FILE, delim = "\t", show_col_types = FALSE,
                    na = c("", "NA", "888", "777", "999"))
trunc <- read_csv(TRUNC_FILE, show_col_types = FALSE, na = c("", "NA")) %>%
  select(-any_of("...1"))

dat <- pheno %>%
  filter(ethnicityClass == 1, !is.na(ageInt)) %>%
  left_join(trunc %>% select(BRIDGES_ID, ends_with("_truncating")), by = "BRIDGES_ID") %>%
  rename_with(~ sub("_truncating", "", .x), ends_with("_truncating")) %>%
  mutate(across(all_of(SIG_GENES), ~ as.integer(replace_na(as.numeric(.x), 0) >= 1)))

controls <- dat %>% filter(status == 0) %>% mutate(outcome = 0L)
dcis     <- dat %>% filter(status == 2, MorphologygroupIndex_corr == "Ductal") %>%
  mutate(outcome = 1L)
d <- bind_rows(controls, dcis) %>% mutate(study = factor(study))

fit_one <- function(gene, use_pl) {
  m <- d %>% select(outcome, carrier = all_of(gene), ageInt, study) %>% drop_na()
  f <- logistf(outcome ~ carrier + ageInt + study, data = m, firth = TRUE, pl = use_pl)
  i <- which(names(coef(f)) == "carrier")
  tibble(gene,
         method  = if (use_pl) "profile-likelihood" else "Wald",
         OR      = exp(coef(f)[i]),
         CI_low  = exp(f$ci.lower[i]),
         CI_high = exp(f$ci.upper[i]))
}

res <- map_dfr(SIG_GENES, function(g) {
  message("  ", g, " — Wald");                w <- fit_one(g, FALSE)
  message("  ", g, " — profile-likelihood");  p <- fit_one(g, TRUE)
  bind_rows(w, p)
})

wide <- res %>%
  mutate(cell = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high)) %>%
  select(gene, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell)

cat("\n=== Wald vs profile-likelihood (DCIS, significant genes) ===\n")
print(as.data.frame(wide))

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write_xlsx(list(comparison = wide, raw = res),
           "outputs/tables/ci_method_check.xlsx")
message("\nSaved: outputs/tables/ci_method_check.xlsx")
