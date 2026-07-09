# =============================================================================
# BMI sensitivity analysis for the Objective-2 hormonal associations
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — University of Cambridge MPhil
#
# PURPOSE
#   Body mass index (BMI) is an established breast-cancer risk factor and a
#   plausible confounder of the hormonal/reproductive associations, but was
#   not in the primary Objective-2 models. This script tests whether the
#   FDR-significant hormonal associations are robust to BMI adjustment. To
#   isolate the effect of adjustment (rather than the change in sample), each
#   association is refitted on the BMI-complete subset both WITHOUT and WITH
#   BMI as a covariate. BMI is also examined as an exposure in its own right.
#
# OUTPUT
#   outputs/tables/bmi_sensitivity.xlsx
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(writexl)})
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                    show_col_types = FALSE, na = c("", "NA", "888", "777", "999"))

eur <- pheno %>% filter(ethnicityClass == 1, !is.na(ageInt))
ctrl <- eur %>% filter(status == 0)
sets <- list(
  BBD  = bind_rows(eur %>% filter(status==0, BBD_history==1) %>% mutate(outcome=1L),
                   eur %>% filter(status==0, BBD_history==0) %>% mutate(outcome=0L)),
  DCIS = bind_rows(ctrl %>% mutate(outcome=0L),
                   eur %>% filter(status==2, MorphologygroupIndex_corr=="Ductal") %>% mutate(outcome=1L)),
  LCIS = bind_rows(ctrl %>% mutate(outcome=0L),
                   eur %>% filter(status==2, MorphologygroupIndex_corr=="Lobular") %>% mutate(outcome=1L))
)

# The FDR-significant Objective-2 associations to stress-test
targets <- tribble(
  ~outcome, ~exposure,   ~orig_OR,
  "BBD",    "parity",    0.94,
  "DCIS",   "parity",    0.91,
  "LCIS",   "parity",    0.75,
  "DCIS",   "menoStat",  0.51,
  "BBD",    "HRTEver",   1.64
)

or_of <- function(d, formula_str, term) {
  fit <- glm(as.formula(formula_str), data = d, family = binomial())
  unname(exp(coef(fit)[term]))
}

res <- pmap_dfr(targets, function(outcome, exposure, orig_OR) {
  d <- sets[[outcome]] %>%
    select(outcome, all_of(exposure), ageInt, study, BMI) %>%
    mutate(study = factor(study)) %>% drop_na()      # BMI-complete subset
  or_no  <- or_of(d, paste("outcome ~", exposure, "+ ageInt + study"), exposure)
  or_bmi <- or_of(d, paste("outcome ~", exposure, "+ ageInt + study + BMI"), exposure)
  tibble(Outcome = outcome, Exposure = exposure,
         `Original OR (full sample)` = orig_OR,
         `n (BMI-complete)` = nrow(d),
         `OR without BMI` = round(or_no, 3),
         `OR with BMI`    = round(or_bmi, 3),
         `Change` = sprintf("%+.1f%%", 100 * (or_bmi - or_no) / or_no))
})

cat("=== BMI sensitivity: hormonal associations, BMI-complete subset ===\n")
print(as.data.frame(res))

# BMI as an exposure in its own right (per 1 unit), adjusted for age + study
bmi_exp <- map_dfr(names(sets), function(o) {
  d <- sets[[o]] %>% select(outcome, BMI, ageInt, study) %>%
    mutate(study=factor(study)) %>% drop_na()
  fit <- glm(outcome ~ BMI + ageInt + study, data = d, family = binomial())
  ci <- suppressMessages(confint.default(fit)["BMI", ])
  tibble(Outcome = o, `n (complete)` = nrow(d),
         `BMI OR (per unit)` = round(exp(coef(fit)["BMI"]), 3),
         `CI low` = round(exp(ci[1]), 3), `CI high` = round(exp(ci[2]), 3),
         `p` = round(summary(fit)$coefficients["BMI","Pr(>|z|)"], 4))
})
cat("\n=== BMI as an exposure (per 1 kg/m^2), adjusted for age + study ===\n")
print(as.data.frame(bmi_exp))

write_xlsx(list(hormonal_BMI_adjusted = res, BMI_as_exposure = bmi_exp),
           "outputs/tables/bmi_sensitivity.xlsx")
message("\nSaved: outputs/tables/bmi_sensitivity.xlsx")
