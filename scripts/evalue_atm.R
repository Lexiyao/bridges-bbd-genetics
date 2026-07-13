# evalue_atm.R
# E-value sensitivity analysis (VanderWeele & Ding, Ann Intern Med 2017) for the
# strongest DCIS association: ATM protein-truncating variants (Table 4.2 / 4.3a).
#
# The E-value is the minimum strength of association, on the risk-ratio scale, that
# an unmeasured confounder (or, here, an unmeasured selection effect) would need to
# have with BOTH the exposure (ATM carriage) and the outcome (DCIS), above and beyond
# the measured covariates (age, contributing study), to fully explain away the
# observed association. DCIS is rare in the source population, so the odds ratio
# approximates the risk ratio and is used directly.
#
# Reproduces: E-value (point) = 7.44, E-value (CI lower bound) = 3.68.

evalue <- function(rr) if (rr <= 1) NA_real_ else rr + sqrt(rr * (rr - 1))

OR <- 3.99; LO <- 2.13; HI <- 7.02   # ATM truncating vs DCIS

cat(sprintf("ATM truncating - DCIS: OR %.2f (95%% CI %.2f-%.2f)\n", OR, LO, HI))
cat(sprintf("  E-value, point estimate     = %.2f\n", evalue(OR)))
cat(sprintf("  E-value, CI bound near null = %.2f\n", evalue(LO)))

# Optional cross-check with the EValue package (exact same result):
if (requireNamespace("EValue", quietly = TRUE)) {
  print(round(EValue::evalues.OR(est = OR, lo = LO, hi = HI, rare = TRUE), 3))
}
