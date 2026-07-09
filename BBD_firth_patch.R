# =============================================================================
# PATCH for run_firth() — add pre-check before Firth model to suppress
# LAPACK error when a gene has 0 carriers in both cases AND controls
# (affects TP53 in subtype analyses)
#
# Replace your existing run_firth() function with this version.
# Everything else in the script stays identical.
# =============================================================================

run_firth <- function(data, gene) {

  n_carrier_case <- sum(data[[gene]] == 1 & data$outcome == 1, na.rm = TRUE)
  n_carrier_ctrl <- sum(data[[gene]] == 1 & data$outcome == 0, na.rm = TRUE)
  n_case         <- sum(data$outcome == 1, na.rm = TRUE)
  n_ctrl         <- sum(data$outcome == 0, na.rm = TRUE)

  # ── Pre-check: skip if no carriers observed in either group ────────────────
  # This prevents Firth from attempting to estimate a coefficient that is
  # entirely unidentified due to absence of the predictor (e.g. TP53 in
  # subtype analyses). Reported as NA in tables and excluded from FDR.
  if ((n_carrier_case + n_carrier_ctrl) == 0) {
    message(sprintf("  %s: 0 carriers in both groups — skipped (reported as NA)", gene))
    return(tibble(gene, n_case, n_ctrl,
                  n_carrier_case = 0L, n_carrier_ctrl = 0L,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "no carriers"))
  }

  model_data <- data %>%
    select(outcome, carrier = all_of(gene), ageInt, study) %>%
    mutate(study = factor(study)) %>%
    drop_na()

  fit <- tryCatch(
    logistf(
      outcome ~ carrier + ageInt + study,
      data  = model_data,
      firth = TRUE,
      pl    = TRUE
    ),
    error = function(e) {
      message(sprintf("  Firth model failed for %s: %s", gene, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) {
    return(tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
                  OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                  p_value = NA_real_, method = "failed"))
  }

  carrier_idx <- which(names(coef(fit)) == "carrier")

  tibble(gene, n_case, n_ctrl, n_carrier_case, n_carrier_ctrl,
         OR      = exp(coef(fit)[carrier_idx]),
         CI_low  = exp(fit$ci.lower[carrier_idx]),
         CI_high = exp(fit$ci.upper[carrier_idx]),
         p_value = fit$prob[carrier_idx],
         method  = "Firth")
}
