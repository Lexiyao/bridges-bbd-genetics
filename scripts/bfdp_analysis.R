# =============================================================================
# Bayesian False-Discovery Probability (BFDP) — DCIS gene associations
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes"
# University of Cambridge MPhil Population Health Sciences — 2025/26
#
# PURPOSE
#   The BRIDGES consortium graded each association's credibility with the
#   Bayesian false-discovery probability of Wakefield (2007). This script
#   reproduces that grading for every DCIS gene-level association reported
#   in Chapter 4, alongside the Benjamini-Hochberg FDR already applied.
#
#   BH controls the false-discovery RATE across a test family; BFDP grades
#   the credibility of each INDIVIDUAL association given a prior on effect
#   size. The two are complementary and are reported together (Table 4.3a).
#
# METHOD (Wakefield 2007, Am J Hum Genet 81:208-227)
#   For an estimated log-OR (theta_hat) with variance V:
#     W    = (log(OR97.5) / 1.96)^2         # prior variance, prior centred at OR=1
#     ABF  = sqrt(V/(V+W)) * exp( theta_hat^2/2 * W/(V+W) )   # approx Bayes factor
#     BFDP = (ABF * pi0) / (ABF * pi0 + (1 - pi0))            # pi0 = prior null prob
#   V is taken from the profile-likelihood CI actually reported for the
#   association: V = ((log(hi) - log(lo)) / (2*1.96))^2. Using the reported
#   CI keeps BFDP internally consistent with the Firth interval in Table 4.2/4.3.
#
#   Primary prior: 97.5th-percentile prior OR = 5, pi0 = 0.5.
#   Sensitivity  : OR97.5 in {3, 5, 8}.
#
# INPUTS  : the five significant DCIS associations plus every other gene,
#           read from the saved result tables (no re-fit needed — BFDP is a
#           deterministic function of the reported OR and CI). Carrier counts
#           and p-values come from outputs/tables/.
# OUTPUTS : outputs/tables/BFDP_analysis_DCIS.csv
#           outputs/tables/BFDP_prior_sensitivity.csv
#           outputs/tables/Table_4_3a_BFDP.csv
#
# DATA GOVERNANCE: outputs are non-disclosive summary statistics only.
# =============================================================================

# -- 0. PACKAGES --------------------------------------------------------------
required <- c("readxl", "writexl", "dplyr")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(readxl); library(writexl); library(dplyr) })

# -- 1. WAKEFIELD BFDP --------------------------------------------------------
# theta_hat = log(OR); V from the reported 95% CI; W from the prior OR97.5.
#   Under H0: theta_hat ~ N(0, V);  under H1: theta_hat ~ N(0, V + W).
#   ABF = f(theta_hat | H0) / f(theta_hat | H1)                 (Wakefield 2007)
#       = sqrt((V+W)/V) * exp( -theta_hat^2/2 * W / (V*(V+W)) )
#   BFDP = P(H0 | data) = (ABF * pi0) / (ABF * pi0 + (1 - pi0)).
bfdp <- function(OR, lo, hi, OR975 = 5, pi0 = 0.5) {
  theta <- log(OR)
  V <- ((log(hi) - log(lo)) / (2 * qnorm(0.975)))^2
  W <- (log(OR975) / qnorm(0.975))^2
  ABF <- sqrt((V + W) / V) * exp(-(theta^2 / 2) * (W / (V * (V + W))))
  (ABF * pi0) / (ABF * pi0 + (1 - pi0))
}

# -- 2. DCIS ASSOCIATIONS -----------------------------------------------------
# Reported DCIS odds ratios and profile-likelihood CIs (Table 4.2 truncating,
# Table 4.3 missense), with carrier counts, p-values and the BH-FDR verdict.
dcis <- tribble(
  ~Gene,    ~Variant,      ~OR,   ~lo,   ~hi,   ~ncar, ~p,       ~BH,
  "ATM",    "truncating",  3.99,  2.13,  7.02,  18,    0.00006,  "BH-sig",
  "CHEK2",  "truncating",  2.66,  1.51,  4.42,  40,    0.00110,  "BH-sig",
  "BRCA2",  "truncating",  2.82,  1.43,  5.17,  14,    0.00410,  "BH-sig",
  "CHEK2",  "missense",    2.01,  1.29,  3.00,  36,    0.00300,  "BH-sig",
  "TP53",   "missense",    3.65,  1.51,  7.76,  9,     0.00600,  "BH-sig",
  "PALB2",  "missense",    2.07,  1.08,  3.64,  NA,    0.03000,  "BH-ns",
  "BARD1",  "missense",    2.23,  1.09,  4.14,  NA,    0.03000,  "BH-ns",
  "BRCA1",  "truncating",  2.49,  0.66,  6.89,  3,     0.16000,  "BH-ns",
  "PALB2",  "truncating",  2.45,  0.73,  6.75,  6,     0.14000,  "BH-ns",
  "TP53",   "truncating",  5.88,  0.55, 39.31,  1,     0.12000,  "BH-ns",
  "RAD51C", "truncating",  2.83,  0.30, 12.34,  1,     0.30000,  "BH-ns",
  "BARD1",  "truncating",  0.73,  0.01,  6.54,  1,     0.85000,  "BH-ns",
  "RAD51D", "truncating",  1.29,  0.01, 10.73,  0,     0.87000,  "BH-ns"
)

dcis <- dcis %>%
  mutate(
    BFDP = round(mapply(bfdp, OR, lo, hi), 3),
    verdict = case_when(
      !is.na(ncar) & ncar <= 1 ~ "not noteworthy (<=1 carrier)",
      BFDP < 0.05              ~ "noteworthy (<0.05)",
      BFDP < 0.10              ~ "borderline (~0.05)",
      TRUE                     ~ "not noteworthy"
    ),
    CI = sprintf("%.2f-%.2f", lo, hi)
  ) %>%
  arrange(BFDP)

# -- 3. PRIOR SENSITIVITY (OR97.5 = 3 / 5 / 8) --------------------------------
sig5 <- dcis %>% filter(BH == "BH-sig")
prior_sens <- sig5 %>%
  transmute(
    Association = paste(Gene, ifelse(Variant == "truncating", "PTV", "missense")),
    OR, CI,
    `BFDP(OR97.5=3)` = round(mapply(function(o,l,h) bfdp(o,l,h,OR975=3), OR, lo, hi), 3),
    `BFDP(OR97.5=5)` = round(mapply(function(o,l,h) bfdp(o,l,h,OR975=5), OR, lo, hi), 3),
    `BFDP(OR97.5=8)` = round(mapply(function(o,l,h) bfdp(o,l,h,OR975=8), OR, lo, hi), 3)
  ) %>%
  mutate(Verdict = ifelse(pmax(`BFDP(OR97.5=3)`, `BFDP(OR97.5=5)`, `BFDP(OR97.5=8)`) <= 0.05,
                          "robust", "prior-sensitive"))

# -- 4. WRITE OUTPUTS ---------------------------------------------------------
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

analysis_out <- dcis %>%
  transmute(Gene, Variant, OR, `95% CI` = sprintf("%.2f-%.2f", lo, hi),
            `Carrier cases` = ncar, p, BFDP, `BFDP verdict` = verdict, `BH/FDR` = BH)
write.csv(analysis_out, "outputs/tables/BFDP_analysis_DCIS.csv", row.names = FALSE)
write.csv(prior_sens,   "outputs/tables/BFDP_prior_sensitivity.csv", row.names = FALSE)

# Table 4.3a: the five BH-significant associations, ordered by BFDP.
table_4_3a <- sig5 %>%
  arrange(BFDP) %>%
  transmute(Gene,
            `Variant class` = ifelse(Variant == "truncating", "Truncating", "Missense"),
            OR, `95% CI` = sprintf("%.2f-%.2f", lo, hi),
            `Carrier cases` = ncar, p, BFDP,
            Interpretation = verdict)
write.csv(table_4_3a, "outputs/tables/Table_4_3a_BFDP.csv", row.names = FALSE)

cat("BFDP analysis complete. Five BH-significant associations:\n")
print(as.data.frame(analysis_out %>% filter(`BH/FDR` == "BH-sig")), row.names = FALSE)
cat("\nPrior sensitivity:\n"); print(as.data.frame(prior_sens), row.names = FALSE)
